# Clusters splicing data.
# XXX not working.

source("git/utils.r")
source("git/unmix/seq/cluster/writeClustersTreeView.r")

output.path = "git/sort_paper/splicing/treeview"

experimentNames = read.tsv("git/unmix/seq/quant/experimentNames.tsv")

r = read.tsv("git/sort_paper/splicing/isoformFractionOfGene.tsv")


descr = paste(r$gene_id, rownames(r))
names(descr) = rownames(r)

r = r[,8:51]


colnames(r) = experimentNames[ colnames(r), "name" ]
r = r[ , order(colnames(r)) ]
r = as.matrix(r)

# Tests for a row not being constant.
non.const.row = function(x) {
  v = var(x, na.rm=TRUE)
  !is.na(v) && (v > 0)
}

r = r[ apply(r, 1, non.const.row) , ]

# XXX
# r1 = matrix(100*runif(6785*44), nrow=6785)
# dimnames(r1) = dimnames(r)
# r.orig = r
# r = r1
# r = r+ r.orig

# do clustering
hr = hcluster(as.matrix(r),
  method="correlation", link="complete", nbproc=7)
hc = hcluster(as.matrix(t(r)),
  method="correlation", link="complete", nbproc=7)

clusters = cutree(hr, k=10)
names(clusters) = rownames(r)
cluster.coloring = rep(hsv(c(1:5) / 6, 1, 0.7), 100)

system(paste("mkdir -p", output.path))
write.clusters.treeview(r, hr, hc,
  clusters, cluster.coloring, 
  paste(output.path, "cluster", sep="/"), descr=rownames(r))


