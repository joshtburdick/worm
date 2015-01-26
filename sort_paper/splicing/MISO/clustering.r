# Clusters MISO splicing data.

source("git/utils.r")
source("git/unmix/seq/cluster/writeClustersTreeView.r")

output.path = "git/sort_paper/splicing/MISO/clustering/"

# one of the clusterings
cl = {
  cl1 = read.table("git/cluster/hierarchical/hier.300.clusters/clusters.tsv",
    sep="\t", header=TRUE, as.is=TRUE)
  cl = cl1[,3]
  names(cl) = cl1[,2]
  cl
}
# median expression of those expression clusters
rr = as.matrix(read.tsv("git/cluster/readRatios.tsv"))
rr = rr[ names(cl), c(1:23) ]
# r.median = aggregate(rr, list(cl), function(x) apply(x, 2, median))
r.median = apply(rr, 2, function(x) by(x, cl, median))
rownames(r.median) = paste("cluster", rownames(r.median))

r = read.tsv("git/sort_paper/splicing/MISO/splicingDiff.tsv")

# event names, restricted to cases nearer the 3' end
event.names = read.tsv("git/sort_paper/splicing/MISO/asEventAnnotation.tsv")
# ??? what cutoff to use here?
event.names = event.names[ event.names$dist.to.3p <= 1000 , ]
r = r[ rownames(r) %in% event.names$id , ]

# Tests for a row not being constant.
non.const.row = function(x) {
  v = var(x, na.rm=TRUE)
  !is.na(v) && (v > 0)
}

# Clusters one type of splicing event.
cluster.splicing = function(r, name) {
  i = apply(r, 1, non.const.row)
  r = r[ i , ]

  descr = paste(event.names[ match(rownames(r), event.names$id), "gene" ],
    event.names[ match(rownames(r), event.names$id), "event.type" ],
    rownames(r))
  descr[ grep("cluster", rownames(r)) ] =
    rownames(r)[ grep("cluster", rownames(r)) ]

  # XXX treating missing values as 0.5, for purposes of clustering
  r1 = r
  r1[ is.na(r1) ] = 0

  # do clustering
  hr = hcluster(as.matrix(r1),
    method="correlation", link="complete", nbproc=7)
  hc = hcluster(as.matrix(t(r1)),
    method="correlation", link="complete", nbproc=7)
  # don't shuffle columns
  hc$order = sort(hc$order)

  clusters = cutree(hr, k=trunc(nrow(r) / 20))
  names(clusters) = rownames(r)
  cluster.coloring = rep(hsv(c(1:5) / 6, 1, 0.7), 200)

  system(paste("mkdir -p", output.path))
  write.clusters.treeview(r, hr, hc,
    clusters, cluster.coloring, 
    paste(output.path, name, sep="/"), descr=descr)

  # add .atr file which doesn't re-order the columns
  write.dummy.atr(paste0(output.path, "/", name, ".atr"),
    ncol(r)) 
}


# putting these all in one clustering
cluster.splicing(rbind(as.matrix(r), r.median), "splicing")


# r1 = r[ r$as.type == "SE" , ]
# z = splicing.table.to.matrix(r1)

