# Generate stats between clusters and clinical variables like TMB and molecular_subtype
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(gplots)
  library(corrplot)
  library(RColorBrewer)
  library(FactoMineR)
  library(factoextra)
  library(gtsummary)
  library(webshot)
})

# parse arguments
option_list <- list(
  make_option(c("--xcell_file"), type = "character", help = "xcell score + cluster file (.tsv)"),
  make_option(c("--who_classification"), type = "character", help = "file with updated WHO classification"),
  make_option(c("--histology"), type = "character", help = "histology file (.tsv)"),
  make_option(c("--tmb_file"), type = "character", help = "TMB file (.tsv)"),
  make_option(c("--output_dir"), type = "character", help = "output directory to write files"),
  make_option(c("--plots_dir"), type = "character", help = "output directory to save plots")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
xcell_file <- opt$xcell_file
who_classification <- opt$who_classification
histology <- opt$histology
tmb_file <- opt$tmb_file

# directories
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = F, recursive = T)
plots_dir <- opt$plots_dir
dir.create(plots_dir, showWarnings = F, recursive = T)

# read in xcell score file
xcell_score <- readr::read_tsv(xcell_file)

# read in histology
histology <- readr::read_tsv(histology)

# read WHO annotation
who_classification <- read_tsv(who_classification)
who_classification <- who_classification %>%
  dplyr::select(Kids_First_Biospecimen_ID, `2021_WHO_Classification`) %>%
  unique()

# add to histology
histology <- histology %>%
  left_join(who_classification) %>%
  filter(sample_id %in% xcell_score$sample_id) %>%
  dplyr::select(any_of(c('Kids_First_Biospecimen_ID', '2021_WHO_Classification', 'composition', 'sample_id', 'experimental_strategy')))

# molecular subtypes > 5 samples
mol_subtypes_to_use <- xcell_score %>%
  dplyr::select(Kids_First_Biospecimen_ID, molecular_subtype) %>%
  unique() %>%
  group_by(molecular_subtype) %>%
  summarise(n = n()) %>%
  filter(n >= 5) %>%
  pull(molecular_subtype)

who_class_to_use <- histology %>%
  dplyr::select(Kids_First_Biospecimen_ID, `2021_WHO_Classification`) %>%
  unique() %>%
  dplyr::group_by(`2021_WHO_Classification`) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n >= 5) %>%
  dplyr::pull(`2021_WHO_Classification`)

# create an annotation file using columns of interest
anno_file <- xcell_score %>%
  filter(molecular_subtype %in% mol_subtypes_to_use) %>%
  dplyr::select(Kids_First_Biospecimen_ID,
                molecular_subtype,
                cluster_assigned) %>%
  ungroup() %>%
  distinct() %>%
  arrange(cluster_assigned, molecular_subtype) %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID")

# 1) chi-square test of independence across xcell clusters vs molecular subtype
capture.output(
  table(anno_file$molecular_subtype, anno_file$cluster_assigned),
  file = file.path(output_dir, "xcell_clusters_vs_subtype_chisquare.txt")
)
capture.output(
  chisq.test(
    table(anno_file$molecular_subtype, anno_file$cluster_assigned)
  ),
  file = file.path(output_dir, "xcell_clusters_vs_subtype_chisquare.txt"),
  append = T
)

# 2) generate balloon plot of molecular subtypes with at least 5 samples in a cluster
dat <- anno_file %>%
  group_by(molecular_subtype, cluster_assigned)  %>%
  summarise(n = n()) %>%
  mutate(nmax = max(n)) %>%
  filter(nmax >= 5) %>%
  ungroup() %>%
  dplyr::select(-c(nmax)) %>%
  spread(key = cluster_assigned, value = n, fill = 0) %>%
  column_to_rownames("molecular_subtype")
pdf(file = file.path(plots_dir, "xcell_clusters_vs_subtype_balloonplot.pdf"))
balloonplot(
  x = as.table(as.matrix(t(dat))),
  main = "xCell clusters vs Molecular subtypes",
  xlab = "",
  ylab = "",
  label = T,
  label.size = 0.5,
  rowmar = 4,
  show.margins = FALSE
)
dev.off()

# 3) generate corrplot of molecular subtypes with at least 5 samples in a cluster
chisq <- chisq.test(dat)
pdf(file = file.path(plots_dir, "xcell_clusters_vs_subtype_corrplot.pdf"))
corrplot(
  chisq$residuals,
  is.cor = FALSE,
  tl.srt = 360,
  tl.offset = 1,
  mar = c(1, 2, 1, 1),
  cl.align.text = "l"
)
dev.off()

# 4) Poisson GLM as alternative independence test

anno_file <- xcell_score %>%
  dplyr::inner_join(histology, 'Kids_First_Biospecimen_ID') %>%
  dplyr::filter(`2021_WHO_Classification` %in% who_class_to_use) %>%
  dplyr::select(Kids_First_Biospecimen_ID,
                `2021_WHO_Classification`,
                cluster_assigned) %>%
  ungroup() %>%
  distinct() %>%
  dplyr::arrange(cluster_assigned, `2021_WHO_Classification`) %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID")

anno_file <- anno_file %>%
  dplyr::mutate(
    cluster_assigned = as.factor(cluster_assigned),
    cluster_assigned = relevel(cluster_assigned, ref = 3),
    `2021_WHO_Classification` = as.factor(`2021_WHO_Classification`),
    `2021_WHO_Classification` = relevel(`2021_WHO_Classification`, ref = 'Pediatric-type diffuse low-grade gliomas, NOS')
  ) %>%
  dplyr::group_by(`2021_WHO_Classification`, cluster_assigned) %>%
  dplyr::mutate(y = n())

modp <- glm(y ~ `2021_WHO_Classification` * cluster_assigned,
            data = anno_file,
            family = poisson)
modp_sum <- summary(modp)

# write summary

sink(file = file.path(output_dir, 'poisson_glm_summary.txt'))
modp_sum
sink()

t1 <- tbl_regression(modp, exponentiate = TRUE)

t1 %>%
  as_gt() %>%
  gt::gtsave(filename = file.path(plots_dir, 'poisson_glm_summary.pdf'))

# 5) Correspondence analysis

anno_file_ca <- anno_file %>%
  dplyr::select(c('2021_WHO_Classification', 'cluster_assigned', 'y')) %>%
  unique(by = '2021_WHO_Classification') %>%
  dplyr::mutate(cluster_assigned = paste0('cluster ', cluster_assigned)) %>%
  tidyr::pivot_wider(names_from = 'cluster_assigned', values_from = 'y') %>%
  tibble::column_to_rownames('2021_WHO_Classification')

anno_file_ca[is.na(anno_file_ca)] <- 0

# write data for reproducibility
write_tsv(anno_file_ca %>% rownames_to_column(" "), file = file.path(output_dir, "coranalysis.tsv"))

pdf(file = file.path(plots_dir, 'coranalysis.pdf'))
res.ca <- FactoMineR::CA(anno_file_ca, ncp = 5, graph = TRUE)
dev.off()

# 6) TMB vs xcell-derived clusters

tmb_file <- readr::read_tsv(tmb_file)

# pull DNA samples for the sample ids in xcell file
histology_DNA <- tmb_file %>%
  inner_join(
    histology %>% dplyr::select(Kids_First_Biospecimen_ID, composition, sample_id),
    by = c("Tumor_Sample_Barcode" = "Kids_First_Biospecimen_ID")
  ) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_DNA" = "Tumor_Sample_Barcode") %>%
  dplyr::select(Kids_First_Biospecimen_ID_DNA, sample_id, composition, tmb)

# pull RNA samples for the same sample ids
histology_RNA <- histology %>%
  filter(sample_id %in% histology_DNA$sample_id) %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  dplyr::select(Kids_First_Biospecimen_ID_RNA, sample_id, composition)

# combine with xcell file
histology_RNA <- histology_RNA %>%
  inner_join(
    xcell_score %>%
      dplyr::select(Kids_First_Biospecimen_ID, cluster_assigned) %>% unique(),
    by = c("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID")
  )

# combine both to map tmb from DNA samples and clusters from RNA samples
tmb_df <- histology_RNA %>%
  inner_join(histology_DNA)

# write data for reproducibility
write_tsv(tmb_df, file = file.path(output_dir, "xcell_clusters_vs_tmb.tsv"))

# 1) output violin plot of TMB per cluster
my_comparisons <- list(c("1", "2"), c("1", "3"), c("2", "3"))
pdf(file.path(plots_dir, "xcell_clusters_vs_tmb.pdf"))
tmb_df$cluster_assigned <- as.factor(tmb_df$cluster_assigned)
p <- ggplot(tmb_df,
            aes(x = cluster_assigned, y = log(tmb + 1), fill = cluster_assigned)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  stat_compare_means(comparison = my_comparisons) +
  stat_compare_means(label.y = 1.5, color = "red") +
  ggpubr::theme_pubr() +
  guides(fill = "none") + xlab("Cluster") + ylab("TMB")
print(p)
dev.off()
