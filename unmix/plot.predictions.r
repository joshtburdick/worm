# Plots some predictions, along with the input data which went into them.

output.dir = "git/unmix/plot/ep.20120724/"

source("R/lineage/tree_utils.r")
source("git/plot/plot_expr.r")

load("git/unmix/image/sort_matrix.Rdata")

# the read depth, with and without correction for sorting purity
source("git/unmix/seq/sortPurityCorrection.r")

# definition of which cell contains which gene
gene.on.off = as.matrix(read.csv("R/unmix/sort_paper/unmix/image/sort.matrix.csv.gz",
  as.is=TRUE, row.names=1))
# compute "negative" portions
gene.on.off = {
  negative.fractions = 1 - gene.on.off
  rownames(negative.fractions) = paste(rownames(gene.on.off), "_minus", sep="")
  rbind(gene.on.off, negative.fractions[-1,])
}

# actual images
load("R/unmix/comp_paper/expr.cell.Rdata")

# the pseudoinverse
source("git/unmix/unmix.r")

# the EP predictions
load("git/unmix/ep.20120724.summary.Rdata")

abplp.reporters = read.table("R/unmix/sort_paper/ABplp.reporters.tsv",
  sep="\t", header=TRUE, as.is=TRUE)
# XXX add in negative reporters
abplp.reporters$color = "blue"
abplp.reporters = rbind(abplp.reporters,
  data.frame(gene=c("ceh-36_minus", "cnd-1_minus", "pha-4_minus", "all"),
    color=rep("red", 4), stringsAsFactors=FALSE))

# some hopefully relatively easy-to-predict genes
enriched.lfe = read.table(
  file="R/unmix/sort_paper/unmix/fraction/enriched.fraction.tsv",
  sep="\t", header=TRUE, row.names=1)

# Plots relative expression of markers in various cells.
# Args:
#   M - cell sorting matrix
#   x.f - the expression of the relevant genes, as enrichment relative
#     to control (adjusted for the number of cells sorted)
#   reporters - the reporters to include (??? maybe not used)
#   root - the root of the lineage tree
#   read.counts - the number of reads for this gene
# Side effects: plots a graph of the sorting fractions, colored appropriately
plot.expr.in.fractions = function(M, x.f, reporters,
  root, read.counts) {

  leaf.cells = c(lin.list[[root]], recursive=TRUE)

  # use weighted average for each leaf cell
  cell.to.leaf = make.cell.to.leaf.matrix(leaf.paths)
  m1 = M[reporters,] %*% cell.to.leaf[,leaf.cells]

  # shading by log-ratio of read counts
  m.shading = log2( (read.counts+1) / (read.counts["all"]+1) ) / 2
  m.shading = m.shading[ reporters ]

  m.shading[ m.shading < -1 ] = -1
  m.shading[ m.shading > 1 ] = 1

  # set up axes
  plot(1,1, xlim=c(1, length(leaf.cells)), ylim=c(0, length(reporters)+1),
    type="n", xaxt="n", yaxt="n", main="", xlab="", ylab="")

  # plot amount of each reporter
  for(i in 1:length(reporters)) {
    gene = reporters[i]
    x1 = 1:length(leaf.cells)
    y.p = m.shading[gene] * m1[gene,]
    y.n = (- m.shading[gene]) * m1[gene,]

    col = if (m.shading[gene] >= 0) rgb(y.p, 0, 0) else rgb(0, 0, y.n)

    segments(x1, i, x1, i+1, lwd=6, lend=2, col=col)
  }

  # names of the reporters
  axis(2, at=1:length(reporters),
    labels=paste(sub("_minus", " -", reporters),
    "(", round(read.counts[reporters]), ")"),
    las=2, cex.axis=0.8)
}

# Plots predictions.
plot.prediction.and.fractions.all = function(reporters, m.cell, x.fraction,
  x.actual, x.prediction, root, read.counts) {
#  root = "ABprp"

#  genes = c("hlh-3", "unc-4", rownames(expr.cell), rownames(enriched.lfe))
  genes = dimnames(ep.summary)[[1]]
  genes = c("acbp-6", "acs-17", "acy-2", "aptf-1", "col-144", "C06G8.1", "C25B8.5",
    "C41G11.1", "ceh-43", "col-105", "cpn-3", "mec-17", "mls-1")

  for(g in sort(unique(intersect(genes, dimnames(ep.summary)[[1]])))) {
    cat(g, "")

    rownames(cell.time.on.off)[rownames(cell.time.on.off)=="NA"] = "P0"
    r = data.frame(cell = as.character(rownames(cell.time.on.off)),
      time.1=cell.time.on.off$on, time.2=cell.time.on.off$off, stringsAsFactors=FALSE)

    # plot actual expression, if it's available
if (g %in% rownames(x.actual)) {
    x1 = x.actual[g,lin.node.names]
    col = rgb(scale.to.unit(x1), 0*x1, 0*x1)
    names(col) = colnames(x.actual)
    r1 = r
    r1$col = col[rownames(cell.time.on.off)]
    plot.segments(r1, main=paste(g, "expression from imaging"), root=root, times=c(-25, 400), lwd=4)
}
else {
  plot.new()
}

  # possibly plot EP prediction
  if (g %in% dimnames(ep.summary)[[1]]) {

    x1 = ep.summary[g,,"per.cell.mean"]

    x.max = max(x1, na.rm=TRUE)       # using the overall maximum for now
    x1[ x1 < 0 ] = 0
    x1[ x1 >= x.max ] = x.max

    col = rgb(scale.to.unit(x1), 0*x1, 0*x1)
    names(col) = colnames(x.prediction)
    r1 = r
    r1$col = col[rownames(cell.time.on.off)]
    plot.segments(r1, main=paste(g, "EP prediction"), root=root, times=c(-25, 400), lwd=4)
  }
  else {
    plot.new()
  }

    # plot prediction
    x1 = x.prediction[g,]

#    x.max = max(x1[ c(lin.list[[root]], recursive=TRUE) ], na.rm=TRUE)
    x.max = max(x1, na.rm=TRUE)         # using the overall maximum for now
    x1[ x1 < 0 ] = 0
    x1[ x1 >= x.max ] = x.max
    col = rgb(scale.to.unit(x1), 0*x1, 0*x1)
    names(col) = colnames(x.prediction)
    r1 = r
    r1$col = col[rownames(cell.time.on.off)]
    plot.segments(r1, main=paste(g, "pseudoinverse prediction"), root=root, times=c(-25, 400), lwd=4)

    # plot the expression going into this
    plot.expr.in.fractions(m.cell, x.fraction[g,], reporters$gene, root, read.counts[g,])
  }
}

plot.it = function() {
  system(paste("mkdir -p", output.dir))

  pdf(paste(output.dir, "/predictions_P0.pdf", sep=""), width=12, height=12.5)
  par(mfrow=c(4,1))
  par(mar=c(0.1,8,3,3)+0.1)
  plot.prediction.and.fractions.all(abplp.reporters, gene.on.off, r,
    expr.cell, x.pseudoinverse, "P0", r)
  dev.off()

  pdf(paste(output.dir, "/predictions_ABplp.pdf", sep=""), width=5, height=12.5)
  par(mfrow=c(4,1))
  par(mar=c(0.1,8,3,3)+0.1)
  plot.prediction.and.fractions.all(abplp.reporters, gene.on.off, r,
    expr.cell, x.pseudoinverse, "ABplp", r)
  dev.off()

  pdf(paste(output.dir, "/predictions_ABprp.pdf", sep=""), width=5, height=12.5)
  par(mfrow=c(4,1))
  par(mar=c(0.1,8,3,3)+0.1)
  plot.prediction.and.fractions.all(abplp.reporters, gene.on.off, r,
    expr.cell, x.pseudoinverse, "ABprp", r)
  dev.off()
}

