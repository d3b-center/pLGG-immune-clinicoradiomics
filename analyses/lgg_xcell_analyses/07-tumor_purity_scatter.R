# scatter plot for tumor purity and immune score

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
})

# parse arguments
option_list <- list(
  make_option(c("--mat"), type = "character",
              help = "input matrix (.rds)"),
  make_option(c("--histology"), type = "character",
              help = "histology file (.tsv)"),
  make_option(c("--plots_dir"), type = "character",
              help = "output directory for plots")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
mat <- opt$mat
histology <- opt$histology
plots_dir <- opt$plots_dir
dir.create(plots_dir, showWarnings = F, recursive = T)

# Load the xCell scores, transpose and subset
xcell_scores <-
  readRDS(mat) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>%
  dplyr::rename(immune_score = `immune score`) %>%
  dplyr::select(Kids_First_Biospecimen_ID, immune_score)


# map cohort participant id to biospecimen from histology
histology_df <-
  file.path(histology) %>%
  read_tsv() %>%
  dplyr::select(
    Kids_First_Biospecimen_ID,
    sample_id,
    tumor_fraction,
    tumor_fraction_RFpurify_ABSOLUTE,
    tumor_fraction_RFpurify_ESTIMATE,
    tumor_fraction_LUMP,
    experimental_strategy
  )

# Join the xcell_scores and histology data frames by 'Kids_First_Biospecimen_ID'
Corr_data <- xcell_scores %>%
  left_join(histology_df, by = "Kids_First_Biospecimen_ID")

selected_data <- Corr_data %>%
  dplyr::select(sample_id, immune_score)

# For Tumor Fraction using WGS only

filtered_data <- histology_df %>%
  filter(
    sample_id %in% selected_data$sample_id &
      experimental_strategy == "WGS" & !is.na(tumor_fraction)
  )


# Compute average scores for multiple specimens per sample_id

joined_data <- selected_data %>%
  left_join(filtered_data,
            by = c("sample_id"),
            suffix = c("", ".filtered"))

joined_data <- joined_data %>% 
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(mean_tumor_fraction = mean(tumor_fraction, na.rm = T),
                   mean_immune_score = mean(immune_score, na.rm = T))

joined_data <- joined_data %>%
  filter(!is.na(mean_tumor_fraction))

# 1. Scatter Plot between immune_score and tumor fraction

p <-
  ggplot(joined_data, aes(x = mean_immune_score, y = mean_tumor_fraction)) +
  geom_point() +
  geom_smooth(method = "lm", fullrange = TRUE) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1))

# Add R-squared coefficient and p-value to the plot
p <-
  p + stat_cor(
    method = "pearson",
    label.x = 0.8,
    label.y = 0.8,
    label.sep = "\n"
  )

# save the file in plot directory
pdf_path <-
  file.path(plots_dir, "immunescore_tumorfraction_correlation_plot.pdf")
ggsave(
  filename = pdf_path,
  plot = p,
  width = 10,
  height = 7
)


# For other three plots subset RFpurity data using Methylation

filtered_data <- histology_df %>%
  filter(sample_id %in% selected_data$sample_id &
           experimental_strategy == "Methylation")

# 2. Plot between immune_score and tumor_fraction_RFpurify_ABSOLUTE

joined_data <- selected_data %>%
  left_join(filtered_data,
            by = "sample_id",
            suffix = c("", ".filtered"))

joined_data <- joined_data %>% 
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(mean_tumor_fraction_RFpurify_ABSOLUTE = mean(tumor_fraction_RFpurify_ABSOLUTE, na.rm = T),
                   mean_immune_score = mean(immune_score, na.rm = T))

joined_data <- joined_data %>%
  filter(!is.na(mean_tumor_fraction_RFpurify_ABSOLUTE) &
           !is.na(experimental_strategy))
# Create Correlation Plot
p <-
  ggplot(joined_data,
         aes(x = immune_score, y = tumor_fraction_RFpurify_ABSOLUTE)) +
  geom_point() +
  geom_smooth(method = "lm", fullrange = TRUE) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1))

# Add R-squared coefficient and p-value to the plot
p <-
  p + stat_cor(
    method = "pearson",
    label.x = 0.8,
    label.y = 0.8,
    label.sep = "\n"
  )

# save the file in plot directory
pdf_path <-
  file.path(plots_dir,
            "immunescore_RFpurify_ABSOLUTE_correlation_plot.pdf")
ggsave(
  filename = pdf_path,
  plot = p,
  width = 10,
  height = 7
)

#3. Plot between immune_score and tumor_fraction_RFpurify_ESTIMATE

joined_data <- selected_data %>%
  left_join(filtered_data,
            by = "sample_id",
            suffix = c("", ".filtered"))

joined_data <- joined_data %>% 
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(mean_tumor_fraction_RFpurify_ESTIMATE = mean(tumor_fraction_RFpurify_ESTIMATE, na.rm = T),
                mean_immune_score = mean(immune_score, na.rm = T))

joined_data <- joined_data %>%
  filter(!is.na(mean_tumor_fraction_RFpurify_ESTIMATE) &
           !is.na(experimental_strategy))

# Create Correlation Plot
p <-
  ggplot(joined_data,
         aes(x = immune_score, y = tumor_fraction_RFpurify_ESTIMATE)) +
  geom_point() +
  geom_smooth(method = "lm", fullrange = TRUE) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1))

# Add R-squared coefficient and p-value to the plot
p <-
  p + stat_cor(
    method = "pearson",
    label.x = 0.8,
    label.y = 0.8,
    label.sep = "\n"
  )

# save the file in plot directory
pdf_path <-
  file.path(plots_dir,
            "immunescore_RFpurify_ESTIMATE_correlation_plot.pdf")
ggsave(
  filename = pdf_path,
  plot = p,
  width = 10,
  height = 7
)

#4. Plot between immune_score and tumor_fraction_LUMP

joined_data <- selected_data %>%
  left_join(filtered_data,
            by = "sample_id",
            suffix = c("", ".filtered"))

joined_data <- joined_data %>% 
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(mean_tumor_fraction_LUMP = mean(tumor_fraction_LUMP, na.rm = T),
                mean_immune_score = mean(immune_score, na.rm = T))

joined_data <- joined_data %>%
  filter(!is.na(tumor_fraction_LUMP) &
           !is.na(experimental_strategy))

# Create Correlation Plot
p <-
  ggplot(joined_data, aes(x = immune_score, y = tumor_fraction_LUMP)) +
  geom_point() +
  geom_smooth(method = "lm", fullrange = TRUE) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1))

# Add R-squared coefficient and p-value to the plot
p <-
  p + stat_cor(
    method = "pearson",
    label.x = 0.8,
    label.y = 0.8,
    label.sep = "\n"
  )

# save the file in plot directory
pdf_path <-
  file.path(plots_dir,
            "immunescore_tumor_fraction_LUMP_correlation_plot.pdf")
ggsave(
  filename = pdf_path,
  plot = p,
  width = 10,
  height = 7
)
