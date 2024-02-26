#!/usr/bin/env Rscript


library(RColorBrewer)
library(R.utils)
source("plotting_funcs.R")


prefix="input.treemix.frq"


for(edge in 1:5){
        print (paste0(prefix, '.', edge))
        plot_tree(cex = 0.8, paste0(prefix, ".", edge))
        pdf(file = paste0(prefix, '.', edge, '.pdf'), width = 10, height = 10)
}

