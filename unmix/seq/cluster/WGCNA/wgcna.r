# Clusters using the WGCNA methodology.

library(WGCNA)

options(stringsAsFactors = FALSE)

read.counts = as.matrix(read.table("git/unmix/seq/quant/readsPerMillion_Murray_050912.tsv",
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

r = t(log2(1 + read.counts))

wnet = blockwiseModules(r, maxBlockSize = 5000,
  power = 6, minModuleSize = 20,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE,
#  colors = hsv(1:255, 1,1),
  saveTOMs = TRUE,
  saveTOMFileBase = "git/unmix/seq/cluster/WGCNA/seq-blockwise",
  verbose = 3)

save(wnet, file="git/unmix/seq/cluster/WGCNA/wnet.Rdata")



