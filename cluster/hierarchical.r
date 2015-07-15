# Hierarchical clustering.

# library("flashClust")
# library("fastcluster")   # ??? is this needed?
library("amap")
library("ctc")

source("git/utils.r")
source("git/unmix/seq/cluster/writeClustersTreeView.r")

system("mkdir -p git/cluster/hierarchical/")

# Reads per million (not including the HS or RNAi data.)
rpm = as.matrix(read.tsv("git/cluster/readsPerMillion.tsv"))
rpm = rpm[ , grep("HS|RNAi", colnames(rpm), invert=TRUE, value=TRUE) ]

# Log-transformed mean and max expression - used as a measure of
# whether we saw something expressed in the FACS data.
mean.expr = apply(log2(rpm+3), 1, mean) - 4

# The read ratios.
r = as.matrix(read.tsv("git/cluster/readRatios.tsv"))

# r = r[c(1:8000),]        # XXX for testing

# Tests for a row not being constant.
non.const.row = function(x) {
  v = var(x, na.rm=TRUE)
  !is.na(v) && (v > 0)
}

r = r[ apply(r, 1, non.const.row) , ]

reads.for.clustering = r  # in case other files shadow this

r.sort.only = r[,c(1:23)]

# maximum enrichment / depletion (just in the sorted fractions)
max.enrichment.depletion = apply(abs(r.sort.only), 1, max)

r.sort.only = r.sort.only[ apply(r.sort.only, 1, non.const.row) , ]

# Treat missing numbers as 0.
# XXX this is arguably not ideal.
r[ is.na(r) ] = 0
r.sort.only[ is.na(r.sort.only) ] = 0

# Get list of genes which have TransgenOme clones
transgenome.genes = {
  tg.table = read.table(
    gzfile("data/genetic/TransgenOmeGenes.tsv.gz"),
    sep="\t", header=TRUE, as.is=TRUE)
  unique(c(tg.table[,5], tg.table[,6]))
}

# Gets short descriptions of genes.
func.descr = {
  fd0 = read.table(gzfile(
    "data/wormbase/c_elegans.PRJNA13758.WS240.functional_descriptions.txt.gz"),
    sep="\t", skip=4, quote="", fill=TRUE, as.is=TRUE)
  fd0 = fd0[,c(2,7)]
  fd0 = fd0[ !duplicated(fd0[,1]) , ]
  fd1 = fd0[,2]
  names(fd1) = fd0[,1]
  fd1[ fd1 == "not known" ] = NA
  fd1
}

# Renumbers a clustering, so that it's in ascending order.
# Args:
#   cl - a clustering, assumed to be consecutive numbers
#     starting with 1
#   ordering - the order the numbers are in
# Returns: that clustering, renumbered to be in order.
renumber.clusters.sorted = function(cl, ordering) {
  gene.names = names(cl)
  cl1 = cl[ ordering ]
  n = max(cl1)

  p = order(order(sapply(1:n, function(i) which(i==cl1)[1])))

  r = p[cl]
  names(r) = gene.names
  r
}


# definition of correlation, which skips missing values
# XXX possibly not used, as it's fairly slow
cor.dist = function(a) {
  a = as.dist( 1 - cor(t(as.matrix(a)), use="na.or.complete") )
#  a[ is.na(a) ] = 2
  a
}

# Does hierarchical clustering.
# Args:
#   r - the dataset
#   r.cluster - subset of data to use for clustering
#   output.name - directory in which to save output
#     (relative to this directory)
#   num.clusters - vector of different number of clusters to include
# Side effects: writes results in that directory,
#   including list of genes and clusters.
# FIXME this should replace .cdt files with symlinks, if possible.
h.cluster = function(r, r.cluster, output.name, num.clusters.list) {

  # do clustering
  hr = hcluster(as.matrix(r.cluster),
    method="correlation", link="complete", nbproc=7)
  hc = hcluster(as.matrix(t(r)),
    method="correlation", link="complete", nbproc=7)

  # XXX slower way
#  hr = hclust(cor.dist(r.cluster))
#  hc = hclust(cor.dist(t(r)))

  # ignore the column ordering
  hc$order = sort(hc$order)

#    paste(output.path.1, "/clusters", sep=""))

  cluster.coloring = rep(hsv(c(1:5) / 6, 1, 0.7), 100)

  for (num.clusters in num.clusters.list) {

    output.path.1 = paste("git/cluster/hierarchical/",
      output.name, ".", num.clusters, ".clusters", sep="")

    system(paste("mkdir -p ", output.path.1))

    # write out clusters, at a threshold to get some number of clusters
    clusters = cutree(hr, k=num.clusters)
print(clusters[1:5])
print(class(clusters[1:5]))
# browser()
    # sort clusters, in numerically increasing order
    # (just for convenience)
    clusters = renumber.clusters.sorted(clusters, hr$order)

    write.table(data.frame(gene=rownames(r.cluster), cluster=clusters),
      file=paste(output.path.1, "clusters.tsv", sep="/"),
      sep="\t", row.names=TRUE, col.names=NA)

    # adding "number of cluster" annotation, whether each
    # gene has a TransgeneOme clone, and the short
    # "functional description" from Wormbase
    fd1 = func.descr[ rownames(r) ]
    fd1[ is.na(fd1) ] = ""

    descr = paste(clusters, rownames(r),
      ifelse(rownames(r) %in% transgenome.genes, "(TG)", ""),
      fd1)

    write.clusters.treeview(as.matrix(r), hr, hc,
      clusters, cluster.coloring, 
      paste(output.path.1, "clusterUnannotated", sep="/"), descr=NULL) # descr)

    # XXX also save the row dendrogram, as it may be
    # useful for plotting
    save(hr, clusters,
      file=paste0(output.path.1, "/dendrogram.Rdata"))

    # use handmade dendrogram for columns
    system(paste("cp git/cluster/dummy_atr_file.tsv ",
      output.path.1, "/cluster.atr", sep=""))
  }
}

# Does clustering with various numbers of clusters.
if (TRUE) {
  n.clusters = c(50,100,150,200,250,300)

  # XXX hacking around issues with cbind...
  me.1 = data.frame("mean expr" = mean.expr[rownames(r.sort.only)],
    "max change" = max.enrichment.depletion[rownames(r.sort.only)],
    check.names=FALSE)
  r.1 = r[rownames(r.sort.only),]
  mean.expr.rso = cbind(me.1, r.1)
  h.cluster(mean.expr.rso, r.sort.only, "hier", n.clusters)

  # ... XXX and gene name differences
  me.1 = data.frame("mean expr" = mean.expr[rownames(r)],
    "max change" = max.enrichment.depletion[rownames(r)],
    check.names=FALSE)
  mean.expr.r = cbind(me.1, r)
  h.cluster(mean.expr.r, r, "hier.ts", n.clusters)
}

