#!/usr/bin/env Rscript

library(argparse)
library(ggplot2)


args = commandArgs(trailingOnly=TRUE)


df <- read.table(args[1], header = F, col.names=c('Chrom', 'Start', 'End', 'Nsites', 'NSNP', 'Heterozygosity'))


# Identify peaks of heterozygosity
#het <- subset(df, df$Nsites >= 600000)
#avg <- mean(het$Heterozygosity)
#std <- sd(het$Heterozygosity)
#threshold <- avg + (2 * std)
#avg
#std
# Define peaks as regions where the heterozygosity is 2 standard deviation above the mean
#peaks <- subset(het, het$Heterozygosity >= threshold)
#output <- sub('\\.txt$', '.peaks.txt', args[1])
#write.table(peaks, file = output, col.names=FALSE, row.names=FALSE, sep ="\t")


# Set color based on chromosome number
df$color <- NA
df$color[df$Chrom == 1 | df$Chrom == 3 | df$Chrom == 5 | df$Chrom == 7 | df$Chrom == 9 | df$Chrom == 11 | df$Chrom == 13 | df$Chrom == 15 | df$Chrom == 17 | df$Chrom == 19 | df$Chrom == 21 | df$Chrom == 23 | df$Chrom == 25 | df$Chrom == 27 | df$Chrom == 29] <- "#0868ac"
df$color[df$Chrom == 2 | df$Chrom == 4 | df$Chrom == 6 | df$Chrom == 8 | df$Chrom == 10 | df$Chrom == 12 | df$Chrom == 14 | df$Chrom == 16 | df$Chrom == 18 | df$Chrom == 20 | df$Chrom == 22 | df$Chrom == 24 | df$Chrom == 26 | df$Chrom == 28] <- "#7a0177"


# Plot genome-wide heterozygosity
het <- subset(df, df$Nsites >= 600000)
p <- ggplot(het, aes(x=Start, y=Heterozygosity))+geom_bar(stat = 'identity', fill = het$color)+facet_grid(~Chrom, scales='free_x', space='free_x', switch = 'x')+theme_classic()+scale_x_continuous(expand=c(0,0))+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+ylab("Heterozygosity")+ylim(0, 0.025)
figure <- sub('\\.txt$', '.pdf', args[1])
ggsave(figure, width = 30, height = 10, units = "cm")


# Plot distribution of heterozygosity
d <- ggplot(het, aes(x=Heterozygosity))+geom_histogram(alpha=.8, fill = '#43a2ca', color = '#0868ac')+theme_bw()+ylab("Count")+xlab("Heterozygosity")+xlim(0, 0.025)
figure <- sub('\\.txt$', '.density.pdf', args[1])
ggsave(figure, width = 10, height = 10, units = "cm")

