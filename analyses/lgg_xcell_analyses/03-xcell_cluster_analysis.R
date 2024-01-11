# Generate stats between clusters and clinical variables like TMB and molecular_subtype
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(gplots)
  library(corrplot)
  library(RColorBrewer)
})

# parse arguments
option_list <- list(
  make_option(c("--xcell_file"), type = "character",
              help = "xcell score + cluster file (.tsv)"),
  make_option(c("--histology"), type = "character",
              help = "histology file (.tsv)"),
  make_option(c("--output_dir"), type = "character",
              help = "output directory to write files"),
  make_option(c("--plots_dir"), type = "character",
              help = "output directory to save plots")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
xcell_file <- opt$xcell_file
histology <- opt$histology
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = F, recursive = T)
plots_dir <- opt$plots_dir
dir.create(plots_dir, showWarnings = F, recursive = T)

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# read in xcell score file
xcell_score <- readr::read_tsv(xcell_file)

# molecular subtypes > 5 samples 
mol_subtypes_to_use <- xcell_score %>%
  dplyr::select(Kids_First_Biospecimen_ID, molecular_subtype) %>%
  unique() %>%
  group_by(molecular_subtype) %>%
  summarise(n = n()) %>%
  filter(n >= 5) %>% 
  pull(molecular_subtype)

# create an annotation file using columns of interest 
anno_file <- xcell_score %>% 
  filter(molecular_subtype %in% mol_subtypes_to_use) %>%
  dplyr::select(Kids_First_Biospecimen_ID, molecular_subtype, cluster_assigned) %>%
  ungroup() %>%
  distinct() %>% 
  arrange(cluster_assigned, molecular_subtype) %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID")

# read in histology file to map RNA biospecimens to DNA biospecimens
histology_df <- readr::read_tsv(opt$histology) %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id)

# read in TMB file
tmb <- readr::read_tsv(file.path(data_dir, "snv-mutation-tmb-coding.tsv")) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Tumor_Sample_Barcode) %>%
  dplyr::select(Kids_First_Biospecimen_ID, tmb) %>% 
  dplyr::left_join(histology_df) %>%
  dplyr::rename(Kids_First_Biospecimen_ID_DNA = Kids_First_Biospecimen_ID)

# combine with other clinical variables  
anno_file <- anno_file %>% 
  rownames_to_column("Kids_First_Biospecimen_ID") %>%
  dplyr::left_join(histology_df) %>%
  dplyr::left_join(tmb)

# 1) output violin plot of TMB per cluster
my_comparisons <- list(c("1", "2"), c("1", "3"), c("2", "3"))
pdf(file.path(plots_dir, "xcell_clusters_vs_tmb.pdf"))
anno_file$cluster_assigned <- as.factor(anno_file$cluster_assigned)
p <- ggplot(anno_file, aes(x = cluster_assigned, y = log(tmb + 1), fill = cluster_assigned)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  stat_compare_means(comparison = my_comparisons) +
  stat_compare_means(label.y = 1.5, color = "red") + 
  ggpubr::theme_pubr() +
  guides(fill = "none") + xlab("Cluster") + ylab("TMB")
print(p)
dev.off()

# 2) chi-square test of independence across xcell clusters vs molecular subtype
capture.output(table(anno_file$molecular_subtype, anno_file$cluster_assigned), file = file.path(output_dir, "xcell_clusters_vs_subtype_chisquare.txt"))
capture.output(chisq.test(table(anno_file$molecular_subtype, anno_file$cluster_assigned)), file = file.path(output_dir, "xcell_clusters_vs_subtype_chisquare.txt"), append = T)

# 3) generate balloon plot of molecular subtypes with at least 5 samples in a cluster
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
balloonplot(x = as.table(as.matrix(t(dat))),
            main = "xCell clusters vs Molecular subtypes", 
            xlab = "", ylab = "",
            label = T, 
            label.size = 0.5, 
            rowmar = 4,
            show.margins = FALSE)
dev.off()

# 4) generate corrplot of molecular subtypes with at least 5 samples in a cluster 
chisq <- chisq.test(dat)
pdf(file = file.path(plots_dir, "xcell_clusters_vs_subtype_corrplot.pdf"))
corrplot(chisq$residuals, is.cor = FALSE, tl.srt = 360, tl.offset = 1, mar = c(1, 2, 1, 1), cl.align.text = "l")
dev.off()
