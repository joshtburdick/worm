# Compares genes enriched in sort fractions with
# genes enriched in Spencer data.

source("git/utils.r")
source("git/data/name_convert.r")
source("git/plot/heatmap.r")
source("git/plot/utils.r")

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

# a = overlap.matrix(bool.matrix.to.list(m.unnormalized > 0),
#   bool.matrix.to.list(spencer.m))[,13:1]
a = tissue.enrich(facs.m, spencer.m)[,13:1]
b = t(gene.set.avg.enrich(spencer.cluster, r[ , rownames(facs.m) ]))[,13:1]


# convert to a table
# XXX slow but simple way of doing this
gec = NULL
for(facs in rownames(a))
  for(tissue in colnames(a)) {
    gec = rbind(gec, data.frame(tissue = tissue,
      facs = facs,
      tissue.enrich = a[facs,tissue],
      gene.enrich = b[facs,tissue],
      stringsAsFactors=FALSE))
  }

pdf("git/sort_paper/FACS/spencerComparison/geneEnrichComparison.pdf",
  width=6, height=5.5)
#  width=11, height=4)
# par(mfrow=c(1,3))
# par(mar=c(9,9.5,3,1)+0.1)

tissues = sort(unique(gec$tissue))
tissue.colors = hsv(c(1:length(tissues))/length(tissues), 0.8, 0.8, 0.8)
names(tissue.colors) = tissues

par(mar=c(5,5,1,1)+0.1)
plot(gec$tissue.enrich, gec$gene.enrich,
  main="",
  xlab=expression("Enrichment of tissue " * (log[2] * " scale")),
  ylab=expression("Enrichment of tissue-specific genes " *
    (log[2] * " scale")),
  col = tissue.colors[ gec$tissue ],
  pch=20)
legend("topleft", legend = tissues, col=tissue.colors,
  cex=0.6, pch=20)

s = cor.test(gec$tissue.enrich, gec$gene.enrich)
# compute p more exactly, using t distribution
p = pt(s$statistic, df=s$parameter, lower.tail=FALSE)

legend("bottomright",
  legend=expr.format(expression("r = " * r * ", " * italic("p") * p),
    list(r = round(s$estimate, 2), p = format.p(p))))

# add regression line
m = lm(gec$gene.enrich ~ gec$tissue.enrich)
abline(m, lwd=2, col="#00000040")

dev.off()

colnames(gec) = c("Tissue", "FACS sample",
  "Tissue enrichment", "Gene enrichment")
write.tsv(gec, file=
  "git/sort_paper/FACS/spencerComparison/geneEnrichComparison.tsv")


