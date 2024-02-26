library(ggplot2)


# Plot admixture

k = 10

samps <- read.table("admixture/input.fam")[,1]
cv_error <- read.table("admixture/cv.error", col.names = c('K', 'Error'))
cv_error

p <- ggplot(cv_error, aes(x=K, y=Error), color  = 'black')+geom_line()+theme_bw()
ggsave('cv_error.pdf', width = 20, height = 20, units = 'cm')

runs <- list()
for (i in 1:k){
runs[[i]] <- read.table(paste0("admixture/input.admixture.", i, ".Q"))
}

pdf('admixture.pdf', width=10, height=10)
par(mfrow=c(3,1))
for (i in 1:k){
barplot(t(as.matrix(runs[[i]])), col=topo.colors(i), ylab="Ancestry", border="black")
}
dev.off()
