#!/usr/bin/env Rscript

library(argparse)
library(ggplot2)


args = commandArgs(trailingOnly=TRUE)

# Mutation rate (change if necessary)
mu <- 0.000000011
# Generation interval (change if necessary)
gen <- 9

msmc <- read.table(args[1], header=T)
figure <- sub('\\.txt$', '.pdf', args[1])

pdf(figure, width = 10, height = 10)
plot(msmc$left_time_boundary/mu*gen, (1/msmc$lambda)/mu, log="x",type="n", xlab="Years ago", ylab="Effective population size", ylim=c(0, 260000))
lines(msmc$left_time_boundary/mu*gen, (1/msmc$lambda)/mu, type = 's', col = 'red')
dev.off()

