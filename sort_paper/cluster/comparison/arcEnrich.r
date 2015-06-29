# Comparison with various "arc"-based datasets, such as
#   the Reece-Hoyes et al 2013 Y1H and PPI data.

source("git/data/name_convert.r")

y1h = read.delim(gzfile("data/tf/y1h/reeceHoyes2013_mmc3.tsv.gz"), as.is=T)
ppi = read.delim(gzfile("data/ppi/reeceHoyes2013_mmc7.tsv.gz"), as.is=T)

# one of the clusterings
cl = {
  cl1 = read.table("git/cluster/hierarchical/hier.300.clusters/clusters.tsv",
    sep="\t", header=TRUE, as.is=TRUE)
  cl = cl1[,3]
  names(cl) = cl1[,2]
  cl
}

# Compares how often ends of an arc are in the same cluster.
# Args:
#   a, b - two lists of gene names
#   cl - the clustering of genes to use
arc.cluster.enrich = function(a, b, cl) {
  r = NULL

  # compute cluster of each gene
  a1 = cl[ rename.gene.name.vector(a) ]
  b1 = cl[ rename.gene.name.vector(b) ]
  i = (!is.na(a1)) & (!is.na(b1))
  a1 = a1[ i ]
  b1 = b1[ i ]
  stopifnot(length(a1) == length(b1))

  n.cl = table(cl)
  background.prob = sum(choose(n.cl, 2)) / (sum(n.cl)^2)
  p.same.cluster = mean(a1==b1)

  c(num.arcs = length(a1),
    background.prob = background.prob,
    p.same.cluster = p.same.cluster,
    times.more.than.random = p.same.cluster / background.prob)
}



y1h.enrich = arc.cluster.enrich(y1h$bait.ORF.name, y1h$prey.ORF.name, cl)
ppi.enrich = arc.cluster.enrich(ppi$bait.ORF.name, ppi$prey.ORF.name, cl)




