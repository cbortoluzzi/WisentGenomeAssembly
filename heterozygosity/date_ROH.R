#!/usr/bin/env Rscript

library(argparse)
library(ggplot2)


args = commandArgs(trailingOnly=TRUE)


# We set a generation interval of 9 years for wisent and American bison and 5 for cattle
generation <- 9
# We set a recombination of 1 cM/Mb (this is of course an approximation)
recombination <- 1


# Read in runs of homozygosity identified by Bortoluzzi et al. 2020
df <- read.table(args[1], header=F, col.names = c("Chromosome", "Start", "End", "WindowsNumber", "Heterozygosity"))

# Retain only ROHs that are least 100 Kb long
df_100Kb <- subset(df, df$WindowsNumber >= 10)
df_100Kb$Length_mb <- (df_100Kb$WindowsNumber * 10000 ) / 1000000


# Calculate age of ROH
df_100Kb$g <- 100 / (2 * df_100Kb$Length_mb)
df_100Kb$Years <- (df_100Kb$g * generation) / recombination
df_100Kb$Generation <- df_100Kb$Years / generation


# Very long ROHs are defined as those segments that originated in the last 200 years
df_100Kb$ROH_class <- 'NA'
df_100Kb$ROH_class[df_100Kb$Years <= 200] <- 'Very_long'

# Long ROHs are defined as those segments that orginated between 200 and 1,000 years ago
df_100Kb$ROH_class[df_100Kb$Years > 200 & df_100Kb$Years <= 1000] <- 'Long'

# Medium ROHs are defined as those segments that originated between 1,000 and 3,000 years ago
df_100Kb$ROH_class[df_100Kb$Years > 1000 & df_100Kb$Years <= 3000] <- 'Medium'

# Short ROHs are defined as those segments that originated more than 3,000 years ago
df_100Kb$ROH_class[df_100Kb$Years > 3000] <- 'Short'


# Save results to file
output <- sub('\\.txt$', '.age.txt', args[1])
write.table(format(df_100Kb, digits=4), output, col.names=T, row.names=F, quote=F, sep="\t")
