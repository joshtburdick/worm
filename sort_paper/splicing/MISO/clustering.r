# Clusters MISO splicing data.

source("git/utils.r")
source("git/unmix/seq/cluster/writeClustersTreeView.r")

output.path = "git/sort_paper/splicing/MISO/clustering/"

r = read.tsv(gzfile("git/unmix/seq/quant/MISO/misoSummary.tsv.gz"))
event.names = read.table("git/unmix/seq/quant/MISO/geneAndAS.tsv",
  quote="", header=TRUE, as.is=TRUE)

# limit sample names slightly
sample.names = sort(grep("^(HS|RNAi)", unique(r$sample),
  value=TRUE, invert=TRUE))
r = r[ r$sample %in% sample.names , ]

# Tests for a row not being constant.
non.const.row = function(x) {
  v = var(x, na.rm=TRUE)
  !is.na(v) && (v > 0)
}

# Converts from tabular form to a matrix.
splicing.table.to.matrix = function(r) {
  events = unique(r$event_name)

  a = matrix(NA, nrow=length(events), ncol=length(sample.names))
  rownames(a) = events
  colnames(a) = sample.names

  # XXX slightly ugly
  i = cbind(match(r$event_name, events),
    match(r$sample, sample.names))

  # center this at zero
  a[i] = 2 * r$miso_posterior_mean - 1

  a
}

# Clusters one type of splicing event.
cluster.splicing = function(splice.data, name) {
  r = splicing.table.to.matrix(splice.data)

  i = apply(r, 1, non.const.row)
  r = r[ i , ]

  descr = paste(event.names[ match(rownames(r), event.names$id), "gene" ],
    rownames(r))

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

for(as.type in unique(r$as.type)) {
  write.status(as.type)
  cluster.splicing(r[ r$as.type==as.type, ], as.type)

}

r1 = r[ r$as.type == "SE" , ]
z = splicing.table.to.matrix(r1)

