# function to create radar plot
suppressPackageStartupMessages({
  library(fmsb)
  library(reshape2)
  library(tidyverse)
})

radar_plot_per_cluster <- function(dat, colors_border, colors_in, title){
  
  # add max and min of each variable to show on the plot
  dat <- rbind(rep(max(dat), ncol(dat)), rep(min(dat), ncol(dat)), dat)
  
  radarchart(df = dat, axistype = 1,
             # custom polygon
             pcol = colors_border, pfcol = colors_in, plwd = 2, plty = 1,
             # custom the grid
             cglcol = "grey", cglty = 1, axislabcol = "grey", cglwd = 0.8,
             # custom labels
             vlcex = 0.8,
             title = title)
}
