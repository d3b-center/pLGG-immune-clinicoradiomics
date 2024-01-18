##### Define bubbleplot function 
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(viridis)
})

bubble_plot <- function(x, topN, title){
  # p <- ggplot(data = x, aes(x = cluster, y = vt_score)) +
  #   geom_point(aes(fill= cell_type), size = 10, alpha = 0.5, colour="black", shape = 21, stroke = 0.5) +
  #   theme_bw() + theme_Publication2() + ggtitle(title)     
  # add colors and shapes
  x <- x %>% 
    arrange(cluster, desc(v_score))
  x <- x %>% 
    group_by(cluster) %>% 
    mutate(id_col = 1:topN,
           colors = viridis_pal(option = "plasma")(5),
           shapes = as.numeric(cluster) + 20)
  x$cell_type_lab <- paste0(x$cell_type,'_',x$cluster)
  p <- ggplot(data = x, aes(x = cluster, y = v_score, shape = cell_type_lab)) +
    theme_bw() + 
    geom_point(aes(fill = cell_type_lab), size = 10, alpha =  0.7, position=position_dodge(width=0.05)) +
    scale_fill_manual(limits = x$cell_type_lab,
                      labels = x$cell_type,
                      values = x$colors) +
    scale_shape_manual(limits = x$cell_type_lab, 
                       labels = x$cell_type,
                       values = x$shapes) +
    labs(shape = "Cell Types", fill = "Cell Types")  +
    ggtitle(title) +
    theme(legend.key.size = unit(1, "cm")) +
    guides(fill = guide_legend(ncol=1, override.aes = list(size=7)))
  
  return(p)
}