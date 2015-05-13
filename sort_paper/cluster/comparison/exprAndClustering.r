# Comparison of clustering with larval expression data
# from Liu et al.
# FIXME: correct the p-values

source("git/utils.r")

# embryonic expression patterns
emb.expr = read.tsv("~/data/image/exprCell.tsv")

# larval expression patterns
larval.expr = read.tsv("data/expression/Larval_expr_Liu_tableS4b.txt.gz")


# Renames the gene names slightly.
rename.genes = function(g) {

  # rename "standard" gene names
  g = sub("^([a-z]+)([0-9][0-9]?.*)$", "\\1-\\2", g)

  # alter some of the clone-based gene names
  g = sub("_", ".", g)

  g
}

rownames(larval.expr) = rename.genes(rownames(larval.expr))

emb.r = cor(t(emb.expr))
larval.r = cor(t(larval.expr))

# Compares correlations between genes in the same cluster, and
# in different clusters.
# Args:
#   cl - a clustering (as a numeric vector indexed by gene name)
#   r - a matrix of correlations
# Returns: list containing 
#   g - the genes that were in common between cl and r
#   r.same, r.different - correlations between genes in the same
#     cluster, or in different clusters
# Side effects: plots histograms of the correlation of genes
#   in the same or different clusters, respectively
cor.and.clustering = function(cl, r, main) {

  # restrict to shared genes
  g = intersect(names(cl), rownames(r))
  cl = cl[ g ]
  r = r[ g, g ]

  # mask out diagonal
  diag(r) = NA

  # divide correlations into "within" or "between" clusters
  s = outer(cl, cl, "==")
  r.same = r[ s ]
  r.different = r[ !s ]
  r.same = r.same[ !is.na(r.same) ]
  r.different = r.different[ !is.na(r.different) ]

  wilcox = wilcox.test(r.same, r.different)

  # plot images
  hist(r.different, breaks=50, col="#c0c0c0",
    xlim=c(-1,1), main=main,
    xlab="Correlation of genes in different clusters")
  hist(r.same, breaks=40, col="#808080",
    xlim=c(-1,1), main=main,
    xlab="Correlation of genes in the same cluster")
  mtext(paste0("Wilcoxon p = ",
    signif(2 * wilcox$p.value, 2)), side=3, cex=0.8)

  list(g = g, r.same = r.same, r.different = r.different,
    wilcox = wilcox)
}

# Plots this comparison for a given clustering.
plot.it = function(clustering.file, r, main) {

  # a clustering (for now, just one of them)
  cl = {
    cl1 = read.table(clustering.file,
      sep="\t", header=TRUE, as.is=TRUE)
    cl = cl1[,3]
    names(cl) = cl1[,2]
    cl
  }

  r = cor.and.clustering(cl, r, main)
}

# cluster.sizes = c(50,100,150,200,250,300)
cluster.sizes = c(300)

pdf("git/sort_paper/cluster/comparison/exprAndClustering.pdf",
  width=7.5, height=6)

par(mfrow=c(2,2))

for(s in cluster.sizes) {
  plot.it(
    paste0("git/cluster/hierarchical/hier.", s, ".clusters/clusters.tsv"),
    emb.r, "Embryonic expression")
}

for(s in cluster.sizes) {
  plot.it(
    paste0("git/cluster/hierarchical/hier.", s, ".clusters/clusters.tsv"),
    larval.r, "Larval expression")
}

dev.off()

