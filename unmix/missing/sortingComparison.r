# Comparison of enrichment from FACS sorting with expression
# measurement from RNA-seq.

# the expression data (from movies)
expr.cell =
  as.matrix(read.table("git/unmix/unmix_comp/data/exprCell.tsv",
    sep="\t", header=TRUE, row.names=1, as.is=TRUE))

# the sort matrix (which is similar, but is thresholded "by hand.")
sortMatrix = as.matrix(read.table("git/unmix/image/sort/sortMatrix.tsv",
  sep="\t", header=TRUE, row.names=1, as.is=TRUE))
sortMatrix =
  sortMatrix[ !(rownames(sortMatrix) %in% c("hlh-16", "irx-1")) , ]

# the expression data, as reads per million
rna.seq = {
  rna.seq.20110922 = as.matrix(read.table(
    "git/unmix/seq/quant/readsPerMillion/readsPerMillion_20110922.tsv",
    sep="\t", header=TRUE, row.names=1, as.is=TRUE))
  rna.seq.20120509 = as.matrix(read.table(
    "git/unmix/seq/quant/readsPerMillion/readsPerMillion_Murray_050912.tsv",
    sep="\t", header=TRUE, row.names=1, as.is=TRUE))
  rna.seq.20120928 = as.matrix(read.table(
    "git/unmix/seq/quant/readsPerMillion/readsPerMillion_092812.tsv",
    sep="\t", header=TRUE, row.names=1, as.is=TRUE))
  cbind(rna.seq.20110922, rna.seq.20120509, rna.seq.20120928)
}

# mapping between rows of the sort matrix, and sequenced fractions
experimentNames = read.table("git/unmix/seq/quant/experimentNames.tsv",
  sep="\t", header=TRUE, as.is=TRUE)

# genes for which we have expression from imaging and sequence
# (note that this only includes 116 / 123 genes)
g = intersect(rownames(expr.cell), rownames(rna.seq))

# experiments (including RNA-seq replicates of cnd-1 and pha-4)
seq.experiments = c(
  "cnd-1 1/4", "cnd-1 8/19", "cnd-1 12/14",
  "pha-4 5/9", "pha-4 9/1", "pha-4 12/9",
  "ceh-26", "ceh-27", "ceh-36", "F21D5.9", "mir-57",
  "mls-2", "pal-1", "ttx-3", "unc-130")

img.experiments = sub(" .*", "", seq.experiments)

# positive and negative expression, from imaging
img.pos = expr.cell[g,] %*% t(sortMatrix[img.experiments,])
img.neg = expr.cell[g,] %*% t(1 - sortMatrix[img.experiments,])

# the same thing, from RNA-seq
seq.pos.names = experimentNames[
  match(paste(seq.experiments, "(+)"), experimentNames$name), "X"]
seq.neg.names = experimentNames[
  match(paste(seq.experiments, "(-)"), experimentNames$name), "X"]
rna.seq.pos = rna.seq[g, seq.pos.names]
rna.seq.neg = rna.seq[g, seq.neg.names]
colnames(rna.seq.pos) = seq.experiments
colnames(rna.seq.neg) = seq.experiments

# adjust for proportion of the embryo included by those
# sequencing experiments
fraction.sorted = apply(sortMatrix, 1, mean)
rna.seq.pos = t(t(rna.seq.pos) * fraction.sorted[img.experiments])
rna.seq.neg = t(t(rna.seq.neg) * (1-fraction.sorted[img.experiments]))

# For each sorted fraction, plots the relative expression
# measured by imaging, compared with from sequencing.
# Side effects: plots the relative expression
plot.imaging.and.seq.expr = function() {

  for(i in 1:length(img.experiments)) {
    gene.name = img.experiments[i]
    seq.name = seq.experiments[i]

    image.r = img.pos[,gene.name] /
      (img.pos[,gene.name] + img.neg[,gene.name])
    seq.r = rna.seq.pos[,seq.name] /
      (rna.seq.pos[,seq.name] + rna.seq.neg[,seq.name])

    plot(image.r, seq.r, xlim=c(0,1), ylim=c(0,1),
      xlab="imaging",
      ylab="RNA-seq",
      main=paste(seq.name, "sorting"),
      pch=20, col="#00000080")
    mtext("relative expression of 116 genes", cex=0.65, line=0.5)
#    mtext("relative expression as pos/(pos+neg)", cex=0.65, line=0.5)
  }
}

# Plots the same thing, but for each gene.
plot.per.gene.comparison = function() {
  for(g1 in g) {
    image.r = img.pos[g1,] /
      (img.pos[g1,] + img.neg[g1,])
    seq.r = rna.seq.pos[g1,] /
      (rna.seq.pos[g1,] + rna.seq.neg[g1,])

    plot(image.r, seq.r, xlim=c(0,1), ylim=c(0,1),
      xlab="imaging",
      ylab="RNA-seq",
      main=g1,
      pch=20, col="#00000080")
#    mtext("relative expression as pos/(pos+neg)", cex=0.65, line=0.5)
  }

}

pdf("git/unmix/missing/sortingComparison.pdf", width=7.5, height=10)
par(mfrow=c(4,3))
par(mar=c(4,4,4,4) + 0.1)
plot.imaging.and.seq.expr()
plot.per.gene.comparison()
dev.off()
