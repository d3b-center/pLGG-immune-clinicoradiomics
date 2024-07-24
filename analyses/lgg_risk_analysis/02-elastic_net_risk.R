suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(plyr)
  library(purrr)
  library(tidyr)
  library(glmnet)
  library(caret)
  library(pROC)
  library(hdi)
  library(ggplot2)
  library(effects)
  library(hydroGOF)
  library(MLmetrics)
})


# parse arguments
option_list <- list(
  make_option(c("--histology_file"), type = "character", help = "histology file for all OpenPedCan samples (.tsv)"),
  make_option(
    c("--risk_file"),
    type = "character",
    default = NULL,
    help = "file with risk scores"
  ),
  make_option(c("--output_dir"), type = "character", help = "output directory"),
  make_option(c("--plots_dir"), type = "character", help = "plots directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
histology_file <- opt$histology_file
risk_file <- opt$risk_file
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = F, recursive = T)
plots_dir <- opt$plots_dir
dir.create(plots_dir, showWarnings = F, recursive = T)

# read histology file
histology_file <- histology_file %>%
  fread()

# read pathways file
lgg_pathways <- readRDS(file.path(output_dir, 'ssgsea_matrix.rds'))

# read imaging risk file and pull corresponding bs identifiers
lgg_clusters <- readxl::read_xlsx(risk_file)
histology_file <- histology_file %>%
  filter(
    cohort_participant_id %in% lgg_clusters$SubjectID,
    experimental_strategy == "RNA-Seq"
  ) %>%
  dplyr::select(
    cohort_participant_id,
    Kids_First_Biospecimen_ID,
    molecular_subtype,
    CNS_region,
    reported_gender,
    race,
    age_at_diagnosis_days,
    RNA_library
  )

hist_risk <- merge(histology_file,
                   lgg_clusters,
                   by.x = 'cohort_participant_id',
                   by.y = 'SubjectID')

lgg_pathways <- lgg_pathways %>%
  t() %>%
  as.data.table(keep.rownames = T) %>%
  dplyr::rename('Kids_First_Biospecimen_ID' = 'rn')

lgg_pathways_risk <- hist_risk %>%
  dplyr::inner_join(lgg_pathways, by = 'Kids_First_Biospecimen_ID') %>%
  dplyr::rename('risk_score' = 'Risk Score', 'risk_group' = 'Risk Group') %>%
  dplyr::mutate(age_at_diagnosis_days = as.numeric(age_at_diagnosis_days)) %>%
  dplyr::mutate(molecular_subtype = gsub('-', '_', molecular_subtype)) %>%
  dplyr::mutate(molecular_subtype = gsub('\\, ', '_', molecular_subtype)) %>%
  dplyr::mutate(molecular_subtype = gsub(' ', '_', molecular_subtype)) %>%
  dplyr::mutate(molecular_subtype = gsub('\\/', '_', molecular_subtype)) %>%
  dplyr::mutate(CNS_region = gsub('-', '_', CNS_region)) %>%
  dplyr::mutate(CNS_region = gsub('\\, ', '_', CNS_region)) %>%
  dplyr::mutate(CNS_region = gsub(' ', '_', CNS_region)) %>%
  dplyr::mutate(CNS_region = gsub('\\/', '_', CNS_region)) %>%
  dplyr::mutate(race = gsub('-', '_', race)) %>%
  dplyr::mutate(race = gsub('\\, ', '_', race)) %>%
  dplyr::mutate(race = gsub(' ', '_', race)) %>%
  dplyr::mutate(race = gsub('\\/', '_', race)) %>%
  dplyr::mutate_if(is.character, as.factor)

## Set training and test samples
set.seed(50)
train_data <- lgg_pathways_risk %>%
  dplyr::filter(Cohort == 'Discovery')
test_data <- lgg_pathways_risk %>%
  dplyr::filter(Cohort == 'Replicate')

pathway_names <- colnames(lgg_pathways)[!names(lgg_pathways) %in% c("Kids_First_Biospecimen_ID")]

covariates <- c(
  'molecular_subtype',
  'CNS_region',
  'reported_gender',
  'race',
  'age_at_diagnosis_days'
)

train_formula <- paste0('risk_score ~ ', paste0(c(covariates, pathway_names), collapse = " + "))


# Dummy code categorical predictor variables
x <- model.matrix(as.formula(train_formula), data = as.data.frame(train_data))
x <- x[, -1]

# Convert the outcome (class) to a numerical variable
y <- train_data$risk_score


### Train and test the model

cv_5 <- caret::trainControl(method = "cv", number = 10)

hit_elnet <- caret::train(
  form = as.formula(train_formula),
  data = as.data.frame(train_data),
  method = "glmnet",
  tuneLength = 100,
  trControl = cv_5
)

get_best_result <- function(caret_fit) {
  best <- which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result <- caret_fit$results[best, ]
  rownames(best_result) <- NULL
  best_result
}

### write out accuracy metrics and tuning plots.
best_result <- get_best_result(hit_elnet)
fwrite(
  best_result,
  file = file.path(output_dir, 'estimated_testing_error_cv.txt'),
  sep = '\t'
)

p <- ggplot(data = hit_elnet) +
  theme(legend.position = "none") +
  theme_classic()
ggsave(
  filename = file.path(plots_dir, 'EL_metrics.pdf'),
  plot = p,
  width = 15
)

### Predict on training data
predictions <- hit_elnet %>% predict(train_data)
### Training model performance metrics
summary <- data.table(
  'Scatter Index' = nrmse(predictions, train_data$risk_score, norm = 'maxmin'),
  'Rsquare' = R2(predictions, train_data$risk_score)
)
fwrite(
  summary,
  file = file.path(output_dir, 'EL_training_prediction_summary.txt'),
  sep = '\t'
)

# Make predictions on the test data
predictions <- hit_elnet %>% predict(test_data)
# Model performance metrics
summary <- data.table(
  'Scatter Index' = nrmse(predictions, test_data$risk_score, norm = 'maxmin'),
  'Rsquare' = R2(predictions, test_data$risk_score)
)

fwrite(
  summary,
  file = file.path(output_dir, 'EL_testing_prediction_summary.txt'),
  sep = '\t'
)

coef_min <- coef(hit_elnet$finalModel, hit_elnet$bestTune$lambda)

result <- data.table('pathways' = colnames(x)[which(coef_min != 0)], 'coefficients' = coef_min@x) %>%
  dplyr::arrange(coefficients)

fwrite(result,
       file = file.path(output_dir, 'full_coefficient_table.txt'),
       sep = '\t')

result_up <- result %>%
  dplyr::arrange(desc(coefficients)) %>%
  dplyr::slice_head(n = 20)

result_down <- result %>%
  dplyr::arrange(coefficients) %>%
  dplyr::slice_head(n = 20)

result_filt <- rbindlist(l = list(result_up, result_down))

# write data for reproducibility
fwrite(result_filt,
       file = file.path(output_dir, 'coef_EL_LGGrisk.tsv'),
       sep = '\t')

p2 <- ggplot(result_filt,
             aes(
               x = reorder(pathways, -coefficients),
               y = coefficients,
               fill = coefficients
             )) +
  geom_bar(stat = "identity") + coord_flip() + theme_bw() +
  xlab("") +
  ylab("GLM coefficient") +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))  +
  scale_x_discrete(
    labels = function(x)
      stringr::str_wrap(x, width = 50)
  ) +
  guides(fill = "none")


ggsave(
  filename = file.path(plots_dir, 'coef_EL_LGGrisk.pdf'),
  plot = p2,
  width = 15,
  device = 'pdf'
)
