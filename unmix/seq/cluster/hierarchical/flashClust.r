# Does hierarchical clustering.

# library("flashClust")
library("amap")
library("ctc")

source("git/utils.r")
source("git/unmix/seq/cluster/writeClustersTreeView.r")

# source("git/unmix/seq/cluster/hierarchical/colorGeneCluster.r")

# Reads in the reads.
# r = read.table("git/unmix/seq/cluster/readsNormalized.tsv",
#   sep="\t", header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE)
r = as.matrix(read.tsv("git/unmix/seq/cluster/readsFACSandTS.tsv"))

# r = r[1:2000,]        # XXX for testing

# remove constant rows
r = r[ apply(r, 1, var) > 0 , ]

r.sort.only = r[,c(1:23)]
r.sort.only = r.sort.only[ apply(r.sort.only, 1, var) > 0 , ]



# Does hierarchical clustering.
# Args:
#   r - the dataset
#   r.cluster - data subsetted (or weighted) for clustering
#   output.path - directory in which to save output
#     (relative to this directory)
#   num.clusters - number of clusters to include
# Side effects: writes results in that directory,
#   including list of genes and clusters.
h.cluster = function(r, r.cluster, output.path, num.clusters=200) {
  output.path.1 = paste("git/unmix/seq/cluster/hierarchical" ,
    output.path, sep="/")
  system(paste("mkdir -p ", output.path.1))

  a = hcluster.treeview(r, r.cluster,
    paste(output.path.1, "/clusters", sep=""))

  # write out clusters, at a threshold to get some number of clusters
  clusters = cutree(a$hr, k=num.clusters)
  write.table(data.frame(gene=rownames(r.cluster), cluster=clusters),
    file=paste(output.path.1, "clusters.tsv", sep="/"),
    sep="\t", row.names=TRUE, col.names=NA)

  # color the gene tree according to those clusters
#  cluster.colors = clusters
#  names(cluster.colors) = rownames(r.cluster)
#  color.gtr(paste(output.path.1, "/clusters", sep=""),
#    cluster.colors)

  # write in TreeView format
  write.treeview(r, r.cluster, clusters,
    rainbow(max(clusters))[sample(max(clusters),max(clusters))],
    paste(output.path.1, "clusters", sep="/"))
  # use handmade dendrogram
  system(paste("cp git/unmix/seq/cluster/readsNormalizedDendrogram.csv ",
    output.path.1, "/clusters.atr", sep=""))
}

# Does hierarchical clustering, and saves TreeView files.
# Args:
#   r - data matrix
#   r.cluster - the data, weighted for clustering
#   base.filename - base filename
#   method - which distance metric
#   link - which clustering method
# Modified from ctc::hclust2treeview() .
# Side effects: writes .atr, .gtr, and .cdt files.
# Returns: row and column clustering
hcluster.treeview = function(r, r.cluster, basefile, method = "correlation", link = "complete") {
  nbproc = 7

  hr <- hcluster(r.cluster, method = method, link = link, nbproc = nbproc)
  # XXX column clustering is arbitrary
  hc <- hcluster(t(r), method = method, link = link, nbproc = nbproc)
  hc$order = sort(hc$order)
#  r2atr(hc, file = paste(basefile, ".atr", sep = ""))
#  r2gtr(hr, file = paste(basefile, ".gtr", sep = ""))
#  r2cdt(hr, hc, r[rownames(r.cluster),],
#    file = paste(basefile, ".cdt", sep = ""))

  list(hr=hr, hc=hc)
}

if (TRUE) {
  # weights to down-weight the timeseries data
#  w1 = rep(1, 61)
#  w1[c(39:61)] = 0.1

  h.cluster(r[rownames(r.sort.only),], r.sort.only, "hier")
  h.cluster(r, r, "hier.ts")
}

# do clustering, and save results (deprecated)
if (FALSE) {
  cluster1 = hcluster.treeview(r,
    "git/unmix/seq/cluster/hierarchical/cluster1",
    method="correlation", link="complete")
  save(cluster1, file="git/unmix/seq/cluster/hierarchical/cluster1.Rdata")

# also, cluster using just the sorted data
cluster2 = hcluster.treeview(as.matrix(r.sort.only),
  "git/unmix/seq/cluster/hierarchical/cluster2",
  method="correlation", link="complete")
save(cluster2, file="git/unmix/seq/cluster/hierarchical/cluster2.Rdata")
}

