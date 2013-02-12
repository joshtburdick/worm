# Clusters using the WGCNA methodology.

library(WGCNA)

source("git/unmix/seq/cluster/WGCNA/hclustMerge.r")

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads=7)

count.path = "git/unmix/seq/quant/readsPerMillion"

r1 = as.matrix(
  read.table(paste(count.path, "readsPerMillion_Murray_050912.tsv", sep="/"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

r2 = as.matrix(
  read.table(paste(count.path, "readsPerMillion_092812.tsv", sep="/"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

r = log2(1 + cbind(r1, r2))
r = r[ , order(colnames(r)) ]
r = r[1:100,]

if (TRUE) {
  wnet = blockwiseModules(t(r), maxBlockSize = 50,
    power = 8, minModuleSize = 10,
    reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE,
#  colors = hsv(1:255, 1,1),
    saveTOMs = TRUE,
    saveTOMFileBase = "git/unmix/seq/cluster/WGCNA/wnet/seq-blockwise",
    verbose = 3)
  save(wnet, file="git/unmix/seq/cluster/WGCNA/wnet/wnetTiny.Rdata")
}

if (FALSE) {
  r = read.table("git/unmix/seq/cluster/readsNormalized.tsv",
    sep="\t", header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE)

# remove constant rows
  r = r[ apply(r, 1, var) > 0 , ]

  system("mkdir -p git/unmix/seq/cluster/WGCNA/wnet.ts/")
  wnet.ts = blockwiseModules(t(r), maxBlockSize = 5000,
    power = 8, minModuleSize = 10,
    reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE,
#  colors = hsv(1:255, 1,1),
    saveTOMs = TRUE,
    saveTOMFileBase = "git/unmix/seq/cluster/WGCNA/wnet.ts/seq-blockwise",
    verbose = 3)
  save(wnet.ts, file="git/unmix/seq/cluster/WGCNA/wnet.ts/wnet.Rdata")
}


# Smoke test of writing out dendrogram files.
write.files = function() {
  load("git/unmix/seq/cluster/WGCNA/wnet/wnetTiny.Rdata")
  x = r

  hr = mergeWGCNAdendrograms(wnet)

  # XXX for now, just using flashClust column clustering
  load("git/unmix/seq/cluster/hierarchical/cluster2.Rdata")
  hc = cluster2$hc

  options(stringsAsFactors = TRUE)
  basefile = "git/unmix/seq/cluster/WGCNA/clusterTest"

  r2atr(hc, file = paste(basefile, ".atr", sep = ""))
  r2gtr(hr, file = paste(basefile, ".gtr", sep = ""), distance="pearson")
  r2cdt(hr, hc, x, file = paste(basefile, ".cdt", sep = ""))

  options(stringsAsFactors = FALSE)
}

