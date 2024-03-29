# Obtain Kaplan-Meier and Cox regression Survival Statistic Results and Plots
# Author: Adam Kraya, adapted from Komal Rathi's and Varun Kesherwani's scripts

suppressPackageStartupMessages({
  library(tidyverse)
  library(survminer)
  library(cowplot)
  library(survival)
  library(ggplot2)
  library(forestmodel)
  library(data.table)
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

# read in tis scores and create high/low groups based on median
tis_scores <-
  file.path(analysis_dir, 'results', 'tis_analysis', 'TIS_scores.tsv') %>%
  fread()
tis_scores <- tis_scores %>%
  dplyr::mutate(tis_score_group = case_when(
    score_avg > median(score_avg) ~ 'High',
    score_avg < median(score_avg) ~ 'Low'
  ))

# map cohort participant id to biospecimen from histology
histology_df <-
  file.path(data_dir, '20230826_release.annotated_histologies_subset.tsv') %>% fread()
histology_df <- histology_df %>%
  dplyr::select(
    Kids_First_Participant_ID,
    cohort_participant_id,
    Kids_First_Biospecimen_ID,
    sample_id,
    OS_days,
    OS_status,
    EFS_days,
    EFS_event_type,
    age_last_update_days,
    age_at_diagnosis_days,
    reported_gender,
    race,
    CNS_region
  )

# join tis information with histologies data
survival_data <- tis_scores %>%
  inner_join(histology_df, by = "Kids_First_Biospecimen_ID")

survival_data <- survival_data %>%
  unique(by = 'Kids_First_Participant_ID')

# code a EFS_status column with the following logic:
survival_data <- survival_data %>%
  mutate(EFS_status = case_when(
    EFS_event_type == 'Not Applicable' ~ 0,
    is.na(EFS_event_type) ~ 0,
    grepl('Recurrence|Progressive|Second|disease', EFS_event_type) ~ 1,
  ))

# Encode OS status
survival_data$OS_status <-
  ifelse(survival_data$OS_status == "DECEASED", 1, 0)

survival_data <- survival_data %>%
  mutate(EFS_days = ifelse(is.na(EFS_days), age_last_update_days, EFS_days))

survival_data <- survival_data %>%
  mutate(OS_days = ifelse(is.na(OS_days), age_last_update_days, OS_days))

survival_data <- survival_data %>%
  mutate(
    OS_days = as.numeric(OS_days),
    OS_status = as.numeric(OS_status),
    EFS_days = as.numeric(EFS_days),
    EFS_status = as.numeric(EFS_status)
  )
########################################## survival analysis OS
# generate the kaplan-meier model
fit <-
  survival::survfit(
    as.formula(
      "survival::Surv(as.numeric(OS_days), OS_status) ~
                               tis_score_group"
    ),
    data = survival_data
  )
pdf(
  file = file.path(plots_dir, "KM_OS_plot_lgg_tis.pdf"),
  height = 8,
  width = 10,
  onefile = FALSE
)
survival_plot <- ggsurvplot(
  fit,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.col = "strata",
  linetype = "strata",
  surv.median.line = "hv",
  ggtheme = theme_minimal(),
  palette = c("#E7B800", "#2E9FDF"),
  legend.title = "TIS Group",
  legend.labs = c("High", "Low")
)

print(survival_plot)
dev.off()

########################################## survival analysis EFS

# generate the kaplan-meier model
fit <-
  survival::survfit(as.formula("survival::Surv(EFS_days, EFS_status) ~ tis_score_group"),
                    data = survival_data)
pdf(
  file = file.path(plots_dir, "KM_EFS_plot_lgg_tis.pdf"),
  height = 8,
  width = 10,
  onefile = FALSE
)
survival_plot <- ggsurvplot(
  fit,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.col = "strata",
  linetype = "strata",
  surv.median.line = "hv",
  ggtheme = theme_minimal(),
  palette = c("#E7B800", "#2E9FDF"),
  legend.title = "TIS Group",
  legend.labs = c("High", "Low")
)
print(survival_plot)
dev.off()

# Fit the Cox model with OS- get a convergence warning due to limited iteration..checked with max iter but same warning message

res_cox <- coxph(Surv(OS_days, OS_status) ~ molecular_subtype +
                   tis_score_group,
                 data = survival_data)

pdf(file = file.path(plots_dir, "KM_OS_forestplot_lgg_tis.pdf"),
    onefile = FALSE)
forest_model(
  res_cox,
  covariates = c('tis_score_group'),
  exponentiate = F,
  limits = log(c(.5, 50))
)
dev.off()

write.table(
  summary(res_cox)$coefficients,
  file = file.path(results_dir, "coxph_OS_model_summary_tis.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# Fit the Cox model with EFS- get a convergence warning due to limited iteration..checked with max iter but same warning message

res_cox <-
  coxph(Surv(EFS_days, EFS_status) ~ molecular_subtype +
          tis_score_group,
        data = survival_data)

pdf(
  file = file.path(plots_dir, "KM_EFS_forestplot_lgg_tis.pdf"),
  height = 8,
  width = 10,
  onefile = FALSE
)
forest_model(
  res_cox,
  covariates = c('tis_score_group'),
  exponentiate = F,
  limits = log(c(.5, 50))
)
dev.off()

write.table(
  summary(res_cox)$coefficients,
  file = file.path(results_dir, "coxph_EFS_model_summary_tis.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

fwrite(survival_data, 
       file.path(results_dir, 'lgat_tis_groups.txt'), 
       sep = '\t')
