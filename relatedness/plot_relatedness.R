#!/usr/bin/env Rscript

library(argparse)
library(ggplot2)


args = commandArgs(trailingOnly=TRUE)


# Read in relatedness report
df <- read.table(args[1], header=T)


# Plot relatedness statistic as matrix
plot <- ggplot(df, aes(x=INDV1, y=INDV2, fill = RELATEDNESS_PHI))+geom_tile()+scale_fill_gradient(low = "white", high = "red") +labs(x = "", y = "", title = "")+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('relatedness_matrix.pdf', width = 20, height = 20, units = 'cm')

