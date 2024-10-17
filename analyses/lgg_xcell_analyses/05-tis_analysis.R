# TIS analysis
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(gplots)
  library(ggplot2)
  library(ggpubr)
  library(preprocessCore)
  library(pheatmap)
  library(RColorBrewer)
})

#### Parse command line options ------------------------------------------------
option_list <- list(
  make_option(c("--count_file"), type = "character", help = "count matrix file (.rds) "),
  make_option(c("--xcell_file"), type = "character", help = "xcell score + cluster file (.tsv)"),
  make_option(c("--who_classification"), type = "character", help = "file with updated WHO classification"),
  make_option(c("--tis_genes"), type = "character", help = "Tumor inflammation signature genes file (.tsv)"),
  make_option(c("--output_dir"), type = "character", help = "output directory to write files"),
  make_option(c("--plots_dir"), type = "character", help = "output directory to save plots")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
count_file <- opt$count_file
xcell_file <- opt$xcell_file
who_classification <- opt$who_classification
tis_genes <- opt$tis_genes
output_dir <- opt$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = F)
plots_dir <- opt$plots_dir
dir.create(plots_dir, recursive = TRUE, showWarnings = F)

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")

# read in required files
tis_genes <- readr::read_tsv(tis_genes)
gene_count <- readRDS(count_file)
sample_cluster <- readr::read_tsv(xcell_file) %>%
  dplyr::select(Kids_First_Biospecimen_ID,
                molecular_subtype,
                cluster_assigned) %>%
  distinct()

# add WHO classification to sample cluster file
who_classification <- read_tsv(who_classification)
who_classification <- who_classification %>%
  dplyr::select(Kids_First_Biospecimen_ID, `2021_WHO_Classification`)
sample_cluster <- sample_cluster %>%
  inner_join(who_classification)

# filter gene count to samples of interest
gene_count <- gene_count %>%
  dplyr::select(sample_cluster$Kids_First_Biospecimen_ID)

# normalize gene count
count.norm <- preprocessCore::normalize.quantiles(x = as.matrix(log2(gene_count + 1)), copy = FALSE)
count.norm <- as.data.frame(count.norm)

# filter by genelist and format
count.norm <- count.norm[rownames(count.norm) %in% tis_genes$Genes, ]
dat_sums <- colSums(count.norm)
dat_means <- colMeans(count.norm)
tis_output <- data.frame(
  subject_id = names(dat_sums),
  score_sum = dat_sums,
  score_avg = dat_means
)
tis_output <- tis_output[order(tis_output$subject_id), ]

# annotate cluster
tis_cluster_annotated <- tis_output %>%
  dplyr::rename(Kids_First_Biospecimen_ID = subject_id) %>%
  dplyr::left_join(sample_cluster)
tis_cluster_annotated$cluster_assigned <- factor(tis_cluster_annotated$cluster_assigned,
                                                 levels = c("1", "2", "3"))

# write output
write_tsv(tis_cluster_annotated, file = file.path(output_dir, "TIS_scores.tsv"))

# summary stats for clusters vs TIS scores violin plot
tis_cluster_annotated %>%
  group_by(cluster_assigned) %>%
  mutate(iqr = boxplot.stats(score_avg)$stats[4] - boxplot.stats(score_avg)$stats[2]) %>%
  summarise(min_value = min(score_avg),
            max_value = max(score_avg),
            lower_whisker = boxplot.stats(score_avg)$stats[2] - 1.5 * iqr,
            upper_whisker = boxplot.stats(score_avg)$stats[4] + 1.5 * iqr,
            median_value = median(score_avg),
            lower_quartile_q1 = quantile(score_avg)[2],
            upper_quartile_q3 = quantile(score_avg)[4]) %>%
  unique() %>%
  write_tsv(file = file.path(output_dir, "xcell_clusters_vs_TIS_stats.tsv"))

# 1) output violin plot of tis average scores
my_comparisons <- list(c("1", "2"), c("1", "3"), c("2", "3"))
pdf(file.path(plots_dir, "TIS_cluster_avg_violinplot.pdf"))
p <- ggplot(
  tis_cluster_annotated,
  aes(x = cluster_assigned, y = score_avg, color = cluster_assigned)
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  stat_compare_means(comparison = my_comparisons) +
  stat_compare_means(label.y = 13, color = "red") +
  ggpubr::theme_pubr() +
  guides(color = "none") + xlab("Cluster") + ylab("TIS score")
print(p)
dev.off()

# output violin plot of tis total scores
# pdf(file.path(plots_dir, "TIS_cluster_sum_violinplot.pdf"))
# p <- ggplot(
#   tis_cluster_annotated,
#   aes(x = cluster_assigned, y = score_sum, color = cluster_assigned)
# ) +
#   geom_violin(trim = FALSE) +
#   geom_boxplot(width = 0.1, fill = "white") +
#   stat_compare_means(comparison = my_comparisons) +
#   stat_compare_means(label.y = 230, color = "red") +
#   ggpubr::theme_pubr() +
#   guides(color = "none") + xlab("Cluster") + ylab("TIS score")
# print(p)
# dev.off()

# 2) output heatmap of tis gene expression (quantile normalized)
# molecular subtypes > 5 samples
mol_subtypes_to_use <- sample_cluster %>%
  dplyr::select(Kids_First_Biospecimen_ID, molecular_subtype) %>%
  unique() %>%
  group_by(molecular_subtype) %>%
  summarise(n = n()) %>%
  filter(n >= 5) %>%
  pull(molecular_subtype)

sample_cluster <- sample_cluster %>%
  dplyr::select(Kids_First_Biospecimen_ID, `2021_WHO_Classification`, cluster_assigned) %>%
  dplyr::arrange(cluster_assigned, `2021_WHO_Classification`) %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% colnames(count.norm))
count.norm <- count.norm %>%
  dplyr::select(sample_cluster$Kids_First_Biospecimen_ID)
sample_cluster$cluster_assigned <- as.character(sample_cluster$cluster_assigned)

# read annotation colors
anno_colors <- rlist::list.load(file = file.path(data_dir, "colors.yaml"))
anno_colors <- lapply(anno_colors, function(x)
  unlist(x))
# anno_colors$molecular_subtype <- anno_colors$molecular_subtype[names(anno_colors$molecular_subtype) %in% sample_cluster$molecular_subtype]
anno_colors$`2021_WHO_Classification` <- anno_colors$`2021_WHO_Classification`[names(anno_colors$`2021_WHO_Classification`) %in% sample_cluster$`2021_WHO_Classification`]

# compute kruskal wallis p-values
identical(colnames(count.norm),
          sample_cluster$Kids_First_Biospecimen_ID)
total <- data.frame()
for (i in 1:nrow(count.norm)) {
  x <- kruskal.test(x = count.norm[i, ] %>% as.numeric(),
                    g = sample_cluster$cluster_assigned)
  df <- data.frame(gene = rownames(count.norm)[i], p_value = x$p.value)
  total <- rbind(total, df)
}
total <- total %>%
  mutate(gene = ifelse(
    test = p_value < 0.01,
    yes = paste0(gene, "**"),
    no = ifelse(p_value < 0.05, yes = paste0(gene, "*"), no = NA)
  ))
rownames(count.norm) <- total$gene

# write out heatmap input for reproducibility
write_tsv(
  count.norm %>% as.data.frame %>% rownames_to_column(" "),
  file = file.path(output_dir, "TIS_heatmap_quantile_norm.tsv")
)

# italicize gene names 
newnames <- lapply(
  rownames(count.norm),
  function(x) bquote(italic(.(x))))

# generate heatmap (quantile normalized)
pheatmap(
  mat = as.matrix(count.norm),
  cellheight = 14,
  fontsize = 10,
  labels_row = as.expression(newnames),
  color = bluered(256),
  annotation_col = sample_cluster %>%
    column_to_rownames("Kids_First_Biospecimen_ID"),
  annotation_colors = anno_colors,
  cluster_cols = F,
  scale = "row",
  show_colnames = F,
  silent = T,
  main = paste0(
    "Tumor Inflammation Signature\n Samples = ",
    ncol(count.norm),
    " | Genes = ",
    nrow(count.norm)
  ),
  filename = file.path(plots_dir, "TIS_heatmap_quantile_norm.pdf"),
  width = 15,
  height = 10
)

