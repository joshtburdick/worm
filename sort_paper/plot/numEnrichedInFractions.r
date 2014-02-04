# Counts of the number of genes enriched in
# particular FACS-sorted samples.

source("git/utils.r")

rpm = read.tsv("git/cluster/readsPerMillion.tsv")

r = read.tsv("git/cluster/readRatios.tsv")

r.sort.only = r[,c(1:23)]

r.sort.only.averaged = cbind(r.sort.only[,c(1:4,8:11,15:21)],
  "cnd-1" = apply(r.sort.only[,c(5:7)], 1, mean),
  "pha-4" = apply(r.sort.only[,c(5:7)], 1, mean),
  "singlets" = apply(r.sort.only[,c(5:7)], 1, mean))

rpm.facs = rpm[ , c(1:23,38:54,65:68) ]

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

# number of fractions a given gene is enriched / depleted in
num.fractions.enriched.depleted.old = {
  n = apply( abs(r.sort.only.averaged) >= log2(3), 1, sum)
  n = n[n>0]
  n
}

num.fractions.enriched.depleted = {
  n = apply( abs(r.sort.only.averaged) >= 2, 1, sum)
  n = n[n>0]
  n
}

cat("number of genes enriched/depleted in at least one fraction =",
  length(num.fractions.enriched.depleted), "\n")
pdf("git/sort_paper/plot/numEnrichedInFractions.pdf")
hist(num.fractions.enriched.depleted,
  main="Number of fractions a gene is enriched/depleted in",
  xlab="number of fractions",
  breaks=max(num.fractions.enriched.depleted)+1, col="grey")
dev.off()

cat("table of number of fractions enriched/depleted:\n")
show(table(num.fractions.enriched.depleted))



