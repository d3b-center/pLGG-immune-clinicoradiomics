# Generate bubble and radar plots of top 5 and bottom 5 immune cell types per cluster
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(reshape2)
  library(ggplot2)
})

# parse arguments
option_list <- list(
  make_option(c("--vtest_scores"), type = "character",
              help = "v-test scores file (.tsv)"),
  make_option(c("--plots_dir"), type = "character",
              help = "output directory to save plots")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
vtest_scores <- opt$vtest_scores
plots_dir <- opt$plots_dir

# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "lgg_xcell_analyses")
dir.create(plots_dir, recursive = TRUE, showWarnings = F)

# source functions
# source(file.path(analysis_dir, "utils", "bubble_plot.R"))
source(file.path(analysis_dir, "utils", "radar_plot_per_cluster.R"))

# read in vtest scores file
vtest_score <- readr::read_tsv(vtest_scores)
n_clusters <- sort(unique(vtest_score$cluster))
colors_border = c(rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9), rgb(0.7,0.5,0.1,0.9))
colors_in = c(rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4), rgb(0.7,0.5,0.1,0.4))

# generate radarplot for vtest scores 

# get top 5 cell types
vtest_score_top5 <- vtest_score %>%
  arrange(cluster, desc(v_score)) %>%
  dplyr::group_by(cluster) %>%
  top_n(n = 5, wt = v_score)

# updated radar plot per cluster 
pdf(file = file.path(plots_dir, "vtest_up_per_cluster.pdf"), onefile = TRUE, width = 8)
for(i in 1:length(n_clusters)){
  dat <- vtest_score_top5 %>%
    filter(cluster == n_clusters[i])
  
  # convert into matrix
  dat_top5 <- acast(dat %>% mutate(cluster = paste0("cluster", cluster)), 
                    cluster~cell_type, value.var = "v_score") %>%
    as.data.frame()
  
  # generate plot
  radar_plot_per_cluster(dat = dat_top5, colors_border = colors_border[i], colors_in = colors_in[i], title = paste0("Cluster", i))
}
dev.off()

# # circular bar plot (upregulated)
# for(i in 1:length(n_clusters)){
#   dat <- vtest_score_top5 %>%
#     filter(cluster == n_clusters[i])
#   pdf(file = file.path(plots_dir, paste0("cluster", n_clusters[i], "_vtest_up.pdf")))
#   p <- ggplot(dat, aes(x = cell_type, y = v_score, fill = cell_type)) +
#     geom_col(position = "dodge") +
#     coord_polar() + 
#     xlab("Cell Type") + ylab("V-test score") +
#     labs(fill = "Cell Type") + 
#     theme_bw() +
#     theme(axis.text.x = element_blank(), 
#           axis.text = element_text(size = 12),
#           legend.position = "right",
#           legend.text = element_text(size = 12),
#           panel.grid.major = element_line(colour = "darkgray"),
#           panel.grid.minor = element_line(colour = "darkgray")) +
#   ggtitle(paste0("Cluster", i))
#   print(p)
#   dev.off()
# }

# get bottom 5 cell types
vtest_score_bottom5 <- vtest_score %>%
  arrange(cluster, desc(v_score)) %>%
  dplyr::group_by(cluster) %>%
  top_n(n = -5, wt = v_score)

# updated radar plot per cluster 
pdf(file = file.path(plots_dir, "vtest_down_per_cluster.pdf"), onefile = TRUE, width = 8)
for(i in 1:length(n_clusters)){
  dat <- vtest_score_bottom5 %>%
    filter(cluster == n_clusters[i])
  
  # convert into matrix
  dat_bottom5 <- acast(dat %>% mutate(cluster = paste0("cluster", cluster)), 
                    cluster~cell_type, value.var = "v_score") %>%
    as.data.frame()
  
  # generate plot
  radar_plot_per_cluster(dat = dat_bottom5, colors_border = colors_border[i], colors_in = colors_in[i], title = paste0("Cluster", i))
}
dev.off()

# # circular bar plot (downregulated)
# for(i in 1:length(n_clusters)){
#   dat <- vtest_score_bottom5 %>%
#     filter(cluster == n_clusters[i])
#   pdf(file = file.path(plots_dir, paste0("cluster", n_clusters[i], "_vtest_down.pdf")))
#   p <- ggplot(dat, aes(x = cell_type, y = v_score, fill = cell_type)) +
#     geom_col(position = "dodge") +
#     coord_polar() + 
#     xlab("Cell Type") + ylab("V-test score") +
#     labs(fill = "Cell Type") + 
#     theme_bw() +
#     theme(axis.text.x = element_blank(), 
#           axis.text = element_text(size = 12),
#           legend.position = "right",
#           legend.text = element_text(size = 12),
#           panel.grid.major = element_line(colour = "darkgray"),
#           panel.grid.minor = element_line(colour = "darkgray")) +
#     ggtitle(paste0("Cluster", i))
#   print(p)
#   dev.off()
# }

# output plot
# pdf(file = file.path(plots_dir, "vtest_bubble_plot.pdf"), onefile = TRUE)
# print(bubble_plot(x = vtest_score_top5, topN = 5, 
#                   title = paste0("Top ", 5, " cell types")))
# print(bubble_plot(x = vtest_score_bottom5, topN = 5, 
#                   title = paste0("Bottom ", 5, " cell types")))
# dev.off()
