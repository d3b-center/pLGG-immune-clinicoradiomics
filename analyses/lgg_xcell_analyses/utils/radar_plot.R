# function to create radar plot
suppressPackageStartupMessages({
  library(fmsb)
  library(reshape2)
  library(tidyverse)
})

radar_plot <- function(dat, output_file){
  
  # add max and min of each variable to show on the plot
  dat <- rbind(rep(max(dat), ncol(dat)), rep(min(dat), ncol(dat)), dat)

  # Color vector
  colors_border = c(rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9), rgb(0.7,0.5,0.1,0.9))
  colors_in = c(rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4), rgb(0.7,0.5,0.1,0.4))
  
  pdf(file = output_file)
  radarchart(df = dat, axistype = 1,
             # custom polygon
             pcol = colors_border, pfcol = colors_in, plwd = 2, plty = 1,
             # custom the grid
             cglcol = "grey", cglty = 1, axislabcol = "grey", cglwd = 0.8,
             # custom labels
             vlcex = 0.8)
  
  # add legend
  legend(x = 1, y = 1, 
         legend = rownames(dat[-c(1,2),]), bty = "n", pch = 20, 
         col = colors_in, text.col = "grey", cex = 1, pt.cex = 3)
  dev.off()
}
