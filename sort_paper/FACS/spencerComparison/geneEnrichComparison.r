# Compares genes enriched in sort fractions with
# genes enriched in Spencer data.

source("git/utils.r")
source("git/data/name_convert.r")
source("git/plot/heatmap.r")

# which cells are in which sort fraction / Spencer experiment
source("git/sort_paper/unmix/sortMatrix.r")
facs.m = m.unnormalized[-15,] > 0
facs.m = facs.m[ sort(rownames(facs.m)) , ]

spencer.m = as.matrix(read.tsv(
  "git/sort_paper/validation/spencerSortMatrix.tsv"))
spencer.m = spencer.m[ sort(rownames(spencer.m)) , ]

# the enrichments
r = as.matrix(read.tsv("git/cluster/readRatios.tsv"))

# average the enrichments for cnd-1 and pha-4
r = cbind(
  "cnd-1"=apply(r[,c("cnd-1 8/19","cnd-1 12/14","cnd-1 1/4")], 1, mean),
  "pha-4"=apply(r[,c("pha-4 5/9","pha-4 9/1","pha-4 12/9")], 1, mean),
  r[ , c(7:38) ])

# information about expression clusters from WormBase
load("git/data/wormbase/expr.cluster.Rdata")
wb.cluster = expr.cluster
wb.cluster$gene = rename.gene.name.vector(wb.cluster$gene)
wb.cluster = unique(wb.cluster)
wb.cluster$group.name = expr.cluster.descr[ wb.cluster$group ]

# just limiting to Spencer "selectively enriched" clusters for now
wb.cluster = wb.cluster[ , c("gene", "group") ]
wb.cluster = wb.cluster[
  grep("WBPaper00037950.*_embryo_SelectivelyEnriched", wb.cluster$group) , ]

# munge names
wb.cluster$group = sub("WBPaper00037950:", "", wb.cluster$group)
wb.cluster$group = sub("_embryo_SelectivelyEnriched", "", wb.cluster$group)

spencer.cluster = list()
for(g in sort(unique(wb.cluster$group)))
  spencer.cluster[[g]] = unique(wb.cluster[wb.cluster$group==g, "gene"])
# XXX match up names
spencer.cluster = spencer.cluster[c(1,3:12,2,13)]
names(spencer.cluster) = rownames(spencer.m)

# Converts a boolean matrix to a list of TRUE elements.
bool.matrix.to.list = function(a) {
  r = NULL
  for(i in rownames(a))
    r[[ i ]] = colnames(a)[ a[i,] ]
  r
}

# Computes average enrichment for various sets of genes.
# Args:
#   gene.set - the gene sets, as a list of character vectors
#   r - expression log ratios
# Returns: matrix of averages of those (for all genes in r)
gene.set.avg.enrich = function(gene.set, r) {
  a = matrix(nrow=length(gene.set), ncol=ncol(r))
  rownames(a) = names(gene.set)
  colnames(a) = colnames(r)

  for(s in names(gene.set)) {
    g = intersect(gene.set[[s]], rownames(r))
    a[s,] = apply(r[g,], 2, mean)
  }

  a
}

# Utility to compute "log-enrichment" of a tissue in
# a set of cells.
le = function(x, C=0.1) {
  log2( (sum(x) + C) / (sum(!x) + C) )
}

# Computes enrichment of Spencer groups of cells (= "tissues")
# in each sort fraction.
tissue.enrich = function(a, b) {
  r = matrix(nrow=nrow(a), ncol=nrow(b))
  rownames(r) = rownames(a)
  colnames(r) = rownames(b)

  for(i in rownames(r))
    for(j in colnames(r))
      r[i,j] = le(b[j, a[i,]]) - le(b[j, !a[i,]])

  r
}

# Jaccard similarity.
# Args: a, b - vectors
# Returns: Jaccard similarity of a and b.
jaccard.sim = function(a, b)
  length(intersect(a,b)) / length(union(a,b))

# Makes a matrix of how much things overlap.
# Args: a, b - two lists of sets of things
# Returns: matrix of proportions overlapping
overlap.matrix = function(a, b) {
  r = matrix(nrow=length(names(a)), ncol=length(names(b)))
  rownames(r) = names(a)
  colnames(r) = names(b)
  for(i in rownames(r))
    for(j in colnames(r))
      r[i,j] = jaccard.sim(a[[i]], b[[j]])
  r
}

pdf("git/sort_paper/FACS/spencerComparison/geneEnrichComparison.pdf",
  width=5.5, height=5.5)
#  width=11, height=4)
# par(mfrow=c(1,3))
# par(mar=c(9,9.5,3,1)+0.1)

# a = overlap.matrix(bool.matrix.to.list(m.unnormalized > 0),
#   bool.matrix.to.list(spencer.m))[,13:1]
a = tissue.enrich(facs.m, spencer.m)[,13:1]
b = t(gene.set.avg.enrich(spencer.cluster, r[ , rownames(facs.m) ]))[,13:1]

# XXX omitting these for now
if (FALSE) {
  plot.heatmap(a, main="Enrichment of tissues",
    cluster=FALSE, z.lim=c(-6, 6),
    col=c(rgb(0,128:0/128,0), rgb(0:128/128,0,0)))
  plot.heatmap(b, main="Enrichment of tissue-specific genes",
    cluster=FALSE, z.lim=c(-2, 2),
    col=c(rgb(0,128:0/128,0), rgb(0:128/128,0,0)))
}

par(mar=c(5,4,4,1)+0.1)
plot(as.vector(a), as.vector(b),
  main="",
  xlab=expression("Enrichment of tissue " * (log[2] * " scale")),
  ylab=expression("Enrichment of tissue-specific genes " *
    (log[2] * " scale")), 
  pch=20, col="#00000080")

a1 = as.vector(a)
b1 = as.vector(b)
s = cor.test(a1, b1)
m = lm(b1 ~ a1)
legend("bottomright",
  legend=paste("r =", round(s$estimate, 2)), bty="n")
# FIXME add p-value?
abline(m, lwd=2, col="#00000040")
dev.off()

