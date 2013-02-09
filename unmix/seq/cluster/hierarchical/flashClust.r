
# library("flashClust")
library("amap")
library("ctc")

r = read.table("git/unmix/seq/cluster/readsNormalized.tsv",
  sep="\t", header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE)

# remove constant rows
r = r[ apply(r, 1, var) > 0 , ]

# r = r[1:500,]     # XXX

# Does hierarchical clustering.
# Args:
#   x - data matrix
#   base.filename - base filename
#   method - which distance metric
#   link - which clustering method
# Modified from ctc::hclust2treeview() .
# Side effects: writes .atr, .gtr, and .cdt files.
# Returns: row and column clustering
hcluster.treeview = function(x, basefile, method = "correlation", link = "complete") {
  nbproc = 7

  hr <- hcluster(x, method = method, link = link, nbproc = nbproc)
  hc <- hcluster(t(x), method = method, link = link, nbproc = nbproc)
  r2atr(hc, file = paste(basefile, ".atr", sep = ""))
  r2gtr(hr, file = paste(basefile, ".gtr", sep = ""))
  r2cdt(hr, hc, x, file = paste(basefile, ".cdt", sep = ""))

  list(hr=hr, hc=hc)
}

cluster1 = hcluster.treeview(r,
  "git/unmix/seq/cluster/hierarchical/cluster1",
  method="correlation", link="complete")
save(cluster1, file="git/unmix/seq/cluster/hierarchical/cluster1.Rdata")

