# Clusters using the WGCNA methodology.

library(WGCNA)

options(stringsAsFactors = FALSE)

count.path = "git/unmix/seq/quant/readsPerMillion"

r1 = as.matrix(
  read.table(paste(count.path, "readsPerMillion_Murray_050912.tsv", sep="/"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

r2 = as.matrix(
  read.table(paste(count.path, "readsPerMillion_092812.tsv", sep="/"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

r = log2(1 + cbind(r1, r2))
r = r[ , order(colnames(r)) ]

if (FALSE) {
  wnet = blockwiseModules(t(r), maxBlockSize = 5000,
    power = 6, minModuleSize = 20,
    reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE,
#  colors = hsv(1:255, 1,1),
    saveTOMs = TRUE,
    saveTOMFileBase = "git/unmix/seq/cluster/WGCNA/wnet/seq-blockwise",
    verbose = 3)
  save(wnet, file="git/unmix/seq/cluster/WGCNA/wnet/wnet.Rdata")
}

