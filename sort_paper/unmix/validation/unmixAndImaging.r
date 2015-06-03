# Comparison of unmixing results and expression data
# from lineage tracing.

source("git/utils.r")

# expression from lineage tracing
expr.cell = read.tsv("~/data/image/exprCell.tsv")

# names of genes used in sorting
sort.genes = colnames(read.tsv("git/cluster/readRatios.tsv"))
sort.genes = intersect(sort.genes, rownames(expr.cell))
sort.genes = sort(c(sort.genes, "cnd-1", "pha-4", "mls-2", "pros-1"))

# the predictions
source("git/sort_paper/unmix/pseudoinverseEnrichment.r")

# (fairly arbitrary) definition of when a gene is on or off
x.on.off = expr.cell >= 2000
rownames(x.on.off) = sub("ceh-26", "pros-1", rownames(x.on.off))
# omit genes which were used in sorting
g = setdiff(rownames(x.on.off), sort.genes)
g = intersect(g, rownames(x.pseudoinverse))
x.on.off = x.on.off[ g, ]

# FIXME use setdiff(g, rownames(x.pseudoinverse)) to check
# if any genes need renaming

# For each gene with imaging data, finds the average enrichment
# of that gene in cells expressing it (according to the imaging
# data), compared to those not expressing it.
# Args:
#   x - the predicted expression of each gene, from unmixing
#   x.on.off - call of whether each gene is expressed in each cell
# Returns: data frame with columns
#   x.on, x.off - average enrichment of gene in cells
#     expressing (or not) each gene
#   x.diff - difference in enrichment
enrich.diff.on.off = function(x.unmix, x.on.off) {
  r = NULL
  for(g in rownames(x.on.off)) {
    x.on = mean(x.unmix[ g, x.on.off[g,] ])
    x.off = mean(x.unmix[ g, !x.on.off[g,] ])
    r = rbind(r, c(x.on = x.on, x.off = x.off, x.diff = x.on - x.off))
  }

  rownames(r) = rownames(x.on.off)
  as.data.frame(r)
}

r = enrich.diff.on.off(x.pseudoinverse, x.on.off)

# save results as table
write.tsv(r, "git/sort_paper/unmix/validation/unmixAndImaging.tsv")

# and plot them
pdf("git/sort_paper/unmix/validation/unmixAndImaging.pdf",
  width=10, height=5)
par(mfrow=c(1,2))

plot(r$x.on, r$x.off, xlim=c(-2,5), ylim=c(-2,5), pch=20, col="#00000040")
abline(0, 1, col="#00000020", lwd=2)

hist(r$x.diff, breaks=50, col="grey")

dev.off()

