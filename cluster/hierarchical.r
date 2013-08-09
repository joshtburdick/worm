# Hierarchical clustering.

# library("flashClust")
# library("fastcluster")   # ??? is this needed?
library("amap")
library("ctc")

source("git/utils.r")
source("git/unmix/seq/cluster/writeClustersTreeView.r")

# The read ratios.
r = as.matrix(read.tsv("git/cluster/readRatios.tsv"))

r = r[1:1000,]        # XXX for testing

# only keep rows which are not constant
non.const.row = function(x) {
  v = var(x, na.rm=TRUE)
  !is.na(v) && (v > 0)
}
r = r[ apply(r, 1, non.const.row) , ]


# was: r = r[ apply(r, 1, var) > 0 , ]
# remove rows which are all NA
# r = r[ apply(!is.na(r), 1, sum) > 0 , ]

r.sort.only = r[,c(1:23)]
# r.sort.only = r.sort.only[ apply(!is.na(r.sort.only), 1, sum) > 0 , ]

r.sort.only = r.sort.only[ apply(r.sort.only, 1, non.const.row) , ]


# Does hierarchical clustering, and saves TreeView files.
# Args:
#   r - data matrix
#   r.cluster - the data, weighted for clustering
#   method - which distance metric
#   link - which clustering method
# Modified from ctc::hclust2treeview() .
# Returns: row and column clustering
hcluster.treeview = function(r, r.cluster, method = "correlation", link = "complete") {
  nbproc = 7

  hr <- hcluster(r.cluster, method = method, link = link, nbproc = nbproc)

  # XXX column clustering is arbitrary
  hc <- hcluster(t(r), method = method, link = link, nbproc = nbproc)
  hc$order = sort(hc$order)

  list(hr=hr, hc=hc)
}

# Does hierarchical clustering.
# Args:
#   r - the dataset
#   r.cluster - data subsetted (or weighted) for clustering
#   output.path - directory in which to save output
#     (relative to this directory)
#   num.clusters - vector of different number of clusters to include
# Side effects: writes results in that directory,
#   including list of genes and clusters.
# Note that this can crash if some cluster is empty.
h.cluster = function(r, r.cluster, output.path, num.clusters.list) {

  # do the clustering
  a = hcluster.treeview(r, r.cluster)
cat("computed a\n")
#    paste(output.path.1, "/clusters", sep=""))

  for (num.clusters in num.clusters.list) {

    output.path.1 = paste("git/cluster/hierarchical", "/",
      output.path, ".", num.clusters, "clusters", sep="")

    system(paste("mkdir -p ", output.path.1))

    # write out clusters, at a threshold to get some number of clusters
    clusters = cutree(a$hr, k=num.clusters)
    write.table(data.frame(gene=rownames(r.cluster), cluster=clusters),
      file=paste(output.path.1, "clusters.tsv", sep="/"),
      sep="\t", row.names=TRUE, col.names=NA)

    write.clusters.treeview(   
      a$hr, a$hc, r.cluster,
      paste(output.path.1, "clusters", sep="/"))

    # write in TreeView format
#    write.clusters.treeview(r, r.cluster, clusters,
#      rainbow(max(clusters))[sample(max(clusters),max(clusters))],
#      paste(output.path.1, "clusters", sep="/"))

    # use handmade dendrogram
#    system(paste("cp git/unmix/seq/cluster/readsNormalizedDendrogram.csv ",
#      output.path.1, "/clusters.atr", sep=""))
    # FIXME this isn't currently working, so we just omit this
    system(paste("rm ", output.path.1, "/clusters.atr", sep=""))
  }
}

# attempt at simpler alternative method
if (TRUE) {
  # do clustering
  hr = hcluster(as.matrix(r),
    method="correlation", link="complete", nbproc=7)
  hc = hcluster(as.matrix(t(r)),
    method="correlation", link="complete", nbproc=7)
  hc$order = sort(hc$order)

  clusters = cutree(hr, k=20)
  names(clusters) = rownames(r)
  cluster.colors = rep(rainbow(12), 100)

  basefile = "git/cluster/temp"
  write.clusters.treeview(as.matrix(r), hr, hc,
    clusters, cluster.colors, basefile)
  # for now, just removing the array clustering, which
  # will cause TreeView to throw an error message, but at
  # least doesn't shuffle the columns.
  # FIXME: copy in a dummy array clustering
  unlink(paste(basefile, ".atr", sep=""))
}



# do clustering, and save results (deprecated)
if (FALSE) {
  cluster1 = hcluster.treeview(r,
    "git/cluster/hierarchical/cluster1",
    method="correlation", link="complete")
  save(cluster1, file="git/cluster/hierarchical/cluster1.Rdata")

  # also, cluster using just the sorted data
  cluster2 = hcluster.treeview(as.matrix(r.sort.only),
    "git/cluster/hierarchical/cluster2",
    method="correlation", link="complete")
  save(cluster2, file="git/cluster/hierarchical/cluster2.Rdata")
}


# for practice
# r1 = r[c(19160:19354),]
# z = hcluster(r1, method="correlation", nbproc=1)


