# Looks for motifs enriched in particular clusters.

library(rpart)   # deprecated

#motifs = read.table(
#  gzfile("git/tf/motif/motifCount/motifs_1kbUpstream.tsv.gz"),
#  sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)

if (FALSE) {

  source("git/tf/motif/motifCount/motif.counts.r")

}

# Looks for enriched motifs, and writes them in a file.
# Args:
#   path - where to find "clusters.tsv" file
#   c - matrix of predictors
#   output.name - what to name the output file
# Side effects: writes a file of results
t.test.clusters = function(path, x, output.name) {
  cl = read.table(paste(path, "/clusters.tsv", sep=""),
    sep="\t", header=TRUE, row.names=1, stringsAsFactors=TRUE)
  rownames(cl) = cl$gene

  g = intersect(rownames(cl), rownames(x))
  x = x[g,]
  cl = cl[g,]

  r = t.test.all.clusters(cl$cluster, x)
  r$p.bh = p.adjust(r$p, method="hochberg")
  r = r[ r$p.bh <= 0.5 & r$t > 0 , ]
#  r = r[ r$p <= 0.01 & r$t > 0 , ]
  r = r[ order(r$p.bh) , ]

  write.table(r,
    file=paste(path, "/", output.name, ".tsv", sep=""),
    sep="\t", row.names=TRUE, col.names=NA)

  r
}

if (FALSE) {

load("git/unmix/seq/cluster/hierarchical/cluster2.Rdata")
a = cutree(cluster2$hr, k=100)  # XXX arbitrary cutoff

clustering = data.frame(gene=names(a), cluster=a)
rownames(clustering) = names(a)
write.table(clustering[order(clustering$cluster),"cluster",drop=FALSE],
  file="git/unmix/seq/cluster/hierarchical/cluster1.tsv",
  sep="\t", col.names=NA)


# only use genes in both tables
g = intersect(rownames(clustering), rownames(motifs))
clustering = clustering[g,]
motifs = motifs[g,]
}

# one sort of regression...
# cat("doing regression\n")
# a = cbind(cluster = factor(clustering$cluster), motifs)
# f = paste("cluster ~ ", paste(colnames(motifs), collapse="+"))
# r = rpart(formula(f), data=a)

# Does many t-tests.
# Args:
#   a, b - two matrices with the same number of columns
# Returns: data frame with columns:
#   name - which column was being tested
#   t, df, p - results from two-sided t-test
t.test.many = function(a, b) {
  r = NULL

  for(j in colnames(a)) {
    m = t.test(a[,j], b[,j])
    r = rbind(r,
      c(name=j, m[[1]], p=m[[3]], m[[2]]))
  }

  # XXX type conversion hack
  data.frame(name = r[,"name"],
    t = as.numeric(r[,"t"]),
    p = as.numeric(r[,"p"]),
    df = as.integer(r[,"df"]))
}

# Computes t-tests for all the clusters.
# Args:
#   cluster - which cluster each gene is in
#   motif - the motifs
# Returns: data frame of results.
t.test.all.clusters = function(cluster, motif) {
  r = NULL

  for(cl in sort(unique(cluster))) {
cat(cl, "")
    i = cluster == cl

  z1 = (motif[i,])
  z2 = (motif[!i,])

    a = t.test.many(motif[i,], motif[!i,])
#    a = t.test.many(motifs[i,], motifs[ sample(which(!i), sum(i)) , ])
    r = rbind(r, cbind(cluster = cl, a))
  }

  r
}

if (FALSE) {
motifClusterEnrichment = {
  r = t.test.all.clusters(clustering$cluster, motifs)
  r$p.bh = p.adjust(r$p, method="hochberg")
  r = r[ r$p.bh <= 0.5 & r$t > 0 , ]
  r = r[ order(r$p.bh) , ]
  r
}

write.table(motifClusterEnrichment, sep="\t", col.names=NA,
  file="git/unmix/seq/cluster/motifClusterEnrichment.tsv")
}

# t.test.clusters("git/unmix/seq/cluster/hierarchical/hier",
#  motif, "motifEnrichment_5kb")
# t.test.clusters("git/unmix/seq/cluster/hierarchical/hier.ts",
#   motif, "motifEnrichment_5kb")

t.test.clusters("git/unmix/seq/cluster/WGCNA/wnet",
  motif, "motifEnrichment_5kb")
t.test.clusters("git/unmix/seq/cluster/WGCNA/wnet.ts",
  motif, "motifEnrichment_5kb")

# z = t.test.many(motif[c(1:2),], motif[c(100:120),])

