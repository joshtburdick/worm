# Counts of the number of genes enriched in
# particular FACS-sorted samples.

library("grid")
library("ggplot2")

source("git/utils.r")

rpm = read.tsv("git/cluster/readsPerMillion.tsv")

r = read.tsv("git/cluster/readRatios.tsv")

r.sort.only = r[,c(1:23)]

r.sort.only.averaged = cbind(r.sort.only[,c(7:21)],
  "cnd-1" = apply(r.sort.only[,c(1:3)], 1, mean),
  "pha-4" = apply(r.sort.only[,c(4:6)], 1, mean),
  "singlets" = apply(r.sort.only[,c(22:23)], 1, mean))

rpm.facs = rpm[ , c(1:21,36:54,65:68) ]

# max. expression of each gene in some sample
max.expr = apply(rpm.facs, 1, max)

cat("num. genes =", length(max.expr), "\n")

# histogram of highest expression of each gene in each sample
pdf("git/sort_paper/plot/numEnrichedInFractions.pdf",
  width=7, height=10)
par(mfrow=c(2,1))
hist(log10(max.expr), breaks=60, col="grey")
abline(v=0, col="#40404080", lwd=2)

cat("number expressed >= 1 in some sample =", sum(max.expr>=1), "\n")

# number of fractions a gene is expressed in
num.fractions = apply(rpm.facs >= 1, 1, sum)

# hist(num.fractions, breaks=max(num.fractions)+1, col="grey")

dev.off()

cat("max. genes in any individual sample =",
  max(apply(rpm.facs>=1, 2, sum)), "\n")

# restrict to genes considered "expressed"
r.sort.only.averaged = r.sort.only.averaged[ max.expr >= 1, ]
r.sort.only.averaged[ is.na(r.sort.only.averaged) ] = 0

num.fractions.enriched.depleted = {
  n = apply( abs(r.sort.only.averaged) >= 2, 1, sum)
  n = n[n>0]
  n
}

cat("number of genes enriched/depleted in at least one fraction =",
  length(num.fractions.enriched.depleted), "\n")
pdf("git/sort_paper/plot/numEnrichedInFractions.pdf",
  width=7.5, height=5)


if (FALSE) {
  # log-scale version of this
  counts = table(num.fractions.enriched.depleted)
  barplot(log10(counts), yaxt="n",
    main="Number of fractions in which a gene\nis enriched or depleted",
    xlab="Number of fractions", ylab="Count",
    col="grey", cex.main=1.5, cex.axis=0.9, cex.lab=1.2)
  axis(2, labels = c(1,10,100,1000), at=c(0,1,2,3))
  abline(h=0)
} else {
# version with subgraph
#hist(num.fractions.enriched.depleted,
#  main="Number of fractions in which a gene\nis enriched or depleted",
#  xlab="Number of fractions",
#  breaks=max(num.fractions.enriched.depleted)+1, col="grey")
#     subplot( hist(rnorm(100)), 15, 1000)


#Any old plot
#a_plot <- ggplot(cars, aes(speed, dist)) + geom_line()

#A viewport taking up a fraction of the plot area
vp <- viewport(width = 0.6, height = 0.6, x = 0.7, y = 0.7)

     set.seed(5689)
     movies <- movies[sample(nrow(movies), 1000), ]
     # Simple examples
#     qplot(rating, data=movies, geom="histogram")
#     qplot(rating, data=movies, weight=votes, geom="histogram")
#  a_plot <- ggplot(movies, aes(x = rating)) + geom_histogram()

a_plot <- ggplot(data.frame(num=num.fractions.enriched.depleted),
  aes(x=num)) + geom_histogram(binwidth=1) + theme_bw()

#Just draw the plot twice
print(a_plot)
print(a_plot + coord_trans(y="sqrt"), vp = vp)


}
dev.off()

cat("table of number of fractions enriched/depleted:\n")
show(table(num.fractions.enriched.depleted))



