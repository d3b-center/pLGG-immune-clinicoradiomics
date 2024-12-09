# Obtain Kaplan-Meier and Cox regression Survival Statistic Results and Plots

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(survminer)
  library(cowplot)
  library(survival)
  library(ggplot2)
  library(forestmodel)
})

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, 'data')
analysis_dir <- file.path(root_dir, "analyses", "lgg_xcell_analyses")
results_dir <-
  file.path(analysis_dir, "results", "survival_analysis")
dir.create(results_dir, showWarnings = F, recursive = T)
plots_dir <- file.path(analysis_dir, "plots", "survival_analysis")
dir.create(plots_dir, showWarnings = F, recursive = T)

# xcell clusters (imaging data has cohort participant id instead of biospecimen id)
xcell_clusters <-
  file.path(analysis_dir,
            "results",
            "xcell_output",
            "xcell_score_cluster.tsv") %>%
  data.table::fread() %>%
  unique(., by = 'Kids_First_Biospecimen_ID')

# map cohort participant id to biospecimen from histology
histology_df <-
  file.path(data_dir,
            "2024-07-10_annotated_histologies_2021_WHO_class.tsv") %>%
  read_tsv()
histology_df <- histology_df %>%
  dplyr::select(
    Kids_First_Participant_ID,
    Kids_First_Biospecimen_ID,
    sample_id,
    OS_days,
    OS_status,
    EFS_days,
    age_last_update_days,
    age_at_diagnosis_days,
    reported_gender,
    race,
    `2021_WHO_Classification`
  )

# Pull in old release for CNS_region
histology_df_cns <-
  file.path(data_dir, "20230826_release.annotated_histologies.tsv") %>%
  read_tsv() %>%
  dplyr::select(any_of(c(
    'Kids_First_Biospecimen_ID', 'CNS_region'
  )))

histology_df <- histology_df %>%
  dplyr::inner_join(histology_df_cns, by = 'Kids_First_Biospecimen_ID') %>%
  dplyr::mutate(
    `2021_WHO_Classification` = ifelse(
      `2021_WHO_Classification` == 'N/A',
      'Pediatric-type diffuse low-grade gliomas, NOS',
      `2021_WHO_Classification`
    )
  )

# join clusters information with histologies data
survival_data <- xcell_clusters %>%
  dplyr::rename("clusterID" = "cluster_assigned") %>%
  dplyr::inner_join(histology_df, by = "Kids_First_Biospecimen_ID")

# code a EFS_status column with the following logic:
survival_data <- survival_data %>%
  dplyr::mutate(
    OS_days = as.numeric(OS_days),
    EFS_days = as.numeric(EFS_days),
    age_at_diagnosis_days = as.numeric(age_at_diagnosis_days)
  )

survival_data <- survival_data %>%
  dplyr::mutate(
    EFS_status = case_when(
      OS_status == "LIVING" & EFS_days == OS_days ~ 1,
      OS_status %in% c("LIVING", "DECEASED") & is.na(EFS_days) ~ 0,
      OS_status %in% c("LIVING", "DECEASED") &
        EFS_days < OS_days ~ 1,
      OS_status == "DECEASED" & EFS_days == OS_days ~ 0,
    )
  )

survival_data <- survival_data %>%
  dplyr::mutate(EFS_days = ifelse(is.na(EFS_days), age_last_update_days, EFS_days)) %>%
  dplyr::mutate(EFS_days = as.numeric(EFS_days))

survival_data <- survival_data %>%
  dplyr::mutate(OS_days = ifelse(is.na(OS_days), age_last_update_days, OS_days)) %>%
  dplyr::mutate(OS_days = as.numeric(OS_days))

# Encode OS status
survival_data$OS_status <-
  ifelse(survival_data$OS_status == "DECEASED", 1, 0)

# Make cluster ID a factor and sort it
survival_data$clusterID <-
  factor(survival_data$clusterID, levels = sort(as.numeric(unique(
    survival_data$clusterID
  ))))

# remove samples with missing values
survival_data_complete <-
  survival_data[complete.cases(survival_data), ] %>%
  unique(by = 'Kids_First_Participant_ID')

########################################## survival analysis OS
# generate the kaplan-meier model
fit <-
  survival::survfit(as.formula("survival::Surv(as.numeric(OS_days), OS_status) ~ clusterID"),
                    data = survival_data_complete)
pdf(
  file = file.path(plots_dir, "KM_OS_plot_lgg_clusters.pdf"),
  height = 8,
  width = 10,
  onefile = FALSE
)
survival_plot <- ggsurvplot(
  fit,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  linetype = "strata",
  surv.median.line = "hv",
  ggtheme = theme_minimal(),
  palette = c("#E7B800", "#2E9FDF", "#8E44AD"),
  legend.title = "ClusterID",
  legend.labs = c("Cluster 1", "Cluster 2", "Cluster 3"),
  xlab = 'Time (Days)',
  ylab = 'Overall Survival Probability'
)

print(survival_plot)
dev.off()

########################################## survival analysis EFS

# generate the kaplan-meier model
fit <-
  survival::survfit(as.formula("survival::Surv(EFS_days, EFS_status) ~ clusterID"),
                    data = survival_data_complete)

# write out data for reproducibility
write_tsv(broom::tidy(fit),
          file = file.path(results_dir, "KM_EFS_plot_lgg_clusters.tsv"))

pdf(
  file = file.path(plots_dir, "KM_EFS_plot_lgg_clusters.pdf"),
  height = 8,
  width = 10,
  onefile = FALSE
)
survival_plot <- ggsurvplot(
  fit,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  linetype = "strata",
  surv.median.line = "hv",
  ggtheme = theme_minimal(),
  palette = c("#E7B800", "#2E9FDF", "#8E44AD"),
  legend.title = "ClusterID",
  legend.labs = c("Cluster 1", "Cluster 2", "Cluster 3"),
  xlab = 'Time (Days)',
  ylab = 'Progression-Free Survival Probability'
)
print(survival_plot)
dev.off()

# Fit the Cox model with OS, remove subtypes with 2 or fewer samples
survival_data_complete <- survival_data_complete %>%
  dplyr::mutate(clusterID = relevel(clusterID, ref = 3))

survival_data_groups <- survival_data_complete %>%
  dplyr::group_by(`2021_WHO_Classification`) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n > 2) %>%
  dplyr::pull(`2021_WHO_Classification`)

survival_data_complete <- survival_data_complete %>%
  dplyr::filter(`2021_WHO_Classification` %in% survival_data_groups) %>%
  dplyr::mutate(`2021_WHO_Classification` = as.factor(`2021_WHO_Classification`)) %>%
  dplyr::mutate(
    `2021_WHO_Classification` = relevel(`2021_WHO_Classification`, ref = 'Pediatric-type diffuse low-grade gliomas, NOS')
  )

survival_data_complete <- survival_data_complete %>%
  dplyr::rename("Age at diagnosis (Days)" = age_at_diagnosis_days, #this variable is a numeric -> label works
                "Biological Sex" = reported_gender, # factor -> label doesn't work!
                "Immune Cluster" = clusterID  # numeric -> label works...
  )

res_cox <- survival::coxph(Surv(OS_days, OS_status) ~ `Age at diagnosis (Days)` +
                             `Biological Sex` +
                             `Immune Cluster`,
                           data = survival_data_complete)

pdf(
  file = file.path(plots_dir, "KM_OS_forestplot_lgg_clusters.pdf"),
  width = 15,
  onefile = FALSE
)
forest_model(
  res_cox,
  covariates = c('`Age at diagnosis (Days)`', '`Biological Sex`', '`Immune Cluster`'),
  exponentiate = T,
  limits = log(c(.5, 50)),
  exclude_infinite_cis = TRUE
)
dev.off()

write.table(
  summary(res_cox)$coefficients,
  file = file.path(results_dir, "coxph_OS_model_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# Fit the Cox model with EFS- get a convergence warning due to limited iteration..checked with max iter but same warning message

res_cox <-
  survival::coxph(
    Surv(EFS_days, EFS_status) ~ `Age at diagnosis (Days)` +
      `Biological Sex` +
      `Immune Cluster` +
      race + 
      CNS_region + 
      `2021_WHO_Classification`,
    data = survival_data_complete
  )

pdf(
  file = file.path(plots_dir, "KM_EFS_forestplot_lgg_clusters.pdf"),
  height = 8,
  width = 13,
  onefile = FALSE
)
forest_model(
  res_cox,
  covariates = c('`Age at diagnosis (Days)`', '`Biological Sex`', '`Immune Cluster`','`2021_WHO_Classification`'),
  exponentiate = T,
  limits = log(c(.5, 50))
)
dev.off()

write.table(
  summary(res_cox)$coefficients,
  file = file.path(results_dir, "coxph_EFS_model_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
