# Clusters using the WGCNA methodology.

library(WGCNA)

source("git/utils.r")
source("git/unmix/seq/cluster/writeClustersTreeView.r")

options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads=7)

r = as.matrix(read.tsv("git/unmix/seq/cluster/readsFACSandTS.tsv"))

# r = r[1:2000,]        # XXX for testing

# remove constant rows
r = r[ apply(r, 1, var) > 0 , ]

r.sort.only = r[,c(1:23)]
r.sort.only = r.sort.only[ apply(r.sort.only, 1, var) > 0 , ]

# old version of this
if (FALSE) {
r = read.table("git/unmix/seq/cluster/readsNormalized.tsv",
    sep="\t", header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE)
# XXX for prototyping
r = r[1:5000,]
# remove constant rows
r = r[ apply(r, 1, var) > 0 , ]
r.sort.only = r[,c(1:21,23,24)]
r.sort.only = r.sort.only[ apply(r.sort.only, 1, var) > 0 , ]
}

# Runs WGCNA clustering on some dataset.
# Args:
#   r - the dataset
#   r.cluster - data subsetted (or weighted) for clustering
#   output.path - directory in which to save output
#     (relative to this directory)
# Side effects: writes results in that directory
run.wgcna = function(r, r.cluster, output.path) {
  output.path = paste("git/unmix/seq/cluster/WGCNA/",
    output.path, sep="/")

  system(paste("mkdir -p", output.path))
  wnet = blockwiseModules(t(r.cluster), maxBlockSize = 5000,
    power = 6, minModuleSize = 10,
    reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE,
#  colors = hsv(1:255, 1,1),
    saveTOMs = TRUE,
    saveTOMFileBase = paste(output.path, "TOM", sep="/"),
    verbose = 3)
  save(wnet, file=paste(output.path, "wnet.Rdata", sep="/"))

  # write out files
  clusters = wnet$colors + 1

  write.table(data.frame(gene=rownames(r.cluster), cluster=clusters),
    file=paste(output.path, "clusters.tsv", sep="/"),
    sep="\t", row.names=TRUE, col.names=NA)

  # write in TreeView format
  write.treeview(r, r.cluster, clusters,
    # first cluster is colored "white", since it's the
    # things which weren't in any cluster
    c("#ffffff", rainbow(max(clusters-1))[sample(max(clusters-1),max(clusters-1))]),
    paste(output.path, "clusters", sep="/"))
  # use handmade dendrogram
  system(paste("cp git/unmix/seq/cluster/readsNormalizedDendrogram.csv ",
    output.path, "/clusters.atr", sep=""))
}


if (TRUE) {

  # first, run on dataset not including the timeseries data
  run.wgcna(r.sort.only, r.sort.only, "wnet")

  # then, run including the timeseries data
  run.wgcna(r, r, "wnet.ts")
}

