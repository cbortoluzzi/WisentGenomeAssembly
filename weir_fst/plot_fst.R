#!/usr/bin/env Rscript

library(argparse)
library(ggplot2)


args = commandArgs(trailingOnly=TRUE)

# Read in eigenvector
df <- read.table(args[1], header = T)

# Plot histogram distribution of variants per 50 Kb window
p <- ggplot(df, aes(x = N_VARIANTS)) + geom_histogram(fill = '#56B4E9', color = '#0072B2')+theme_bw()+xlab("Number of bi-allelic SNPs per 50 Kb")+ggtitle(args[2])
p + geom_vline(aes(xintercept = 100), color="black", linetype="dashed", size=1)

name_p <- sub('\\.windowed.weir.fst$', '', args[1])
figure <- paste0(name_p, ".pdf")
ggsave(figure, width = 15, height = 15, units = "cm")



