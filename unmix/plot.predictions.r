# Plots some predictions, along with the input data which went into them.

load("R/unmix/sort_paper/unmix/run/ep.20120621.Rdata")
ep.prediction = m

source("R/lineage/tree_utils.r")
source("R/plot/plot_expr.r")

source("R/unmix/sort_paper/unmix/run/pseudoinverse.20120615.r")
x.p[ is.na(x.p) ] = 0
x.p[ x.p < 0 ] = 0
load("R/unmix/sort_paper/unmix/run/pseudoinverse.20120626.Rdata")
trunc.pseudo.p[ trunc.pseudo.p < 0 ] = 0

output.dir = "R/unmix/sort_paper/unmix/plot/ep.and.pseudoinverse/"

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
#   reporters - the reporters to plot, and their colors
#   m - cell sorting matrix
#   x.f - the expression of the relevant genes, as enrichment relative
#     to control (adjusted for the number of cells sorted)
#   root - the root of the lineage tree
#   read.counts - the number of reads for this gene
# Side effects: plots a graph of the sorting fractions, colored appropriately
plot.expr.in.fractions = function(m, x.f, reporters,
  root, read.counts) {
# print(x.f)
  leaf.cells = c(lin.list[[root]], recursive=TRUE)

# m1 = m[reporters$gene,leaf.cells]
# print(dim(m1))
  # use weighted average for each leaf cell
  cell.to.leaf = make.cell.to.leaf.matrix(leaf.paths)
  m1 = m[reporters$gene,] %*% cell.to.leaf[,leaf.cells]

  # shading by log-ratio of read counts
  m.shading = log2( (read.counts+1) / (read.counts["all"]+1) ) / 2
  m.shading = m.shading[ reporters$gene ]

  m.shading[ m.shading < -1 ] = -1
  m.shading[ m.shading > 1 ] = 1
# print(m.shading)
  # set up axes
  plot(1,1, xlim=c(1, length(leaf.cells)), ylim=c(0, length(reporters$gene)+1),
    type="n", xaxt="n", yaxt="n", main="", xlab="", ylab="")
  for(i in 1:nrow(reporters)) {
    gene = reporters[i,"gene"]
    x1 = 1:length(leaf.cells)
#    y1 = m.shading[i] * m1[reporters[i,"gene"],] +
#      (1 - m.shading[i]) * (1 - m1[reporters[i,"gene"]])
    y.p = m.shading[gene] * m1[gene,]
    y.n = (- m.shading[gene]) * m1[gene,]
#    y.p[ y.p < 0.5 ] = 0
#    y.n[ y.n < 0.5 ] = 0
#cat(reporters[i,"gene"], ": ")
#cat("m.shading[i] = ", m.shading[gene], "\n")
#cat("max y.p = ", max(y.p), "\n")
#cat("\n")
    col = if (m.shading[gene] >= 0) rgb(y.p, 0, 0) else rgb(0, 0, y.n)

    segments(x1, i, x1, i+1, lwd=6, lend=2, col=col)   # rgb(y.p, 0, y.n))

#      col=rgb(y1 %*% t(col2rgb(reporters[i,"color"])) / 255))
#     col=rgb(m.shading[i] * m1[reporters[i,"gene"],] %*% t(col2rgb(reporters[i,"color"])) / 255))
#      col=rgb(m1[reporters[i,"gene"],] %*% t(col2rgb(reporters[i,"color"])) / 255))
  }

  # names of the reporters
  axis(2, at=1:nrow(reporters),
    labels=paste(sub("_minus", " -", reporters$gene),
    "(", round(read.counts[reporters$gene]), ")"),
    las=2, cex.axis=0.8)
}

# Plots predictions.
plot.prediction.and.fractions.all = function(reporters, m.cell, x.fraction,
  x.actual, x.prediction, root, read.counts) {
#  root = "ABprp"

  # various neurotransmitters and receptors, etc.
  genes = c("ceh-10","ttx-3","ceh-23","sra-11","kal-1","hen-1","unc-17", "ser-2",
    "ser-1", "ser-4", "ser-5", "ser-7", "dop-1", "dop-2", "dop-3", "dop-4", "dop-5",
    "dop-6", "octr-1", "ser-3", "ser-6", "ser-2", "tyra-2", "tyra-3", "lgc-55",
    "F59D12.1", "T21B4.4", "sprr-1", "sprr-2", "sprr-3",
    "gar-1", "gar-2", "gar-3", "lev-1", "unc-29", "acr-2", "acr-3")

  genes = c("hlh-3", "unc-4", rownames(expr.cell), rownames(enriched.lfe))
#  genes = rownames(ep.prediction)

#  for(g in intersect(rownames(x.actual), rownames(x.prediction))) {
#  for(g in intersect(colnames(x.fraction), intersect(rownames(x.actual), rownames(x.prediction)))) {
  for(g in sort(unique(intersect(genes, rownames(x.prediction))))) {
    cat(g, "")

    rownames(cell.time.on.off)[rownames(cell.time.on.off)=="NA"] = "P0"
    r = data.frame(cell = as.character(rownames(cell.time.on.off)),
      time.1=cell.time.on.off$on, time.2=cell.time.on.off$off, stringsAsFactors=FALSE)

    # plot actual expression, if it's available
if (g %in% rownames(x.actual)) {
    x1 = x.actual[g,lin.node.names]
#    x.max = max(x1[ c(lin.list[[root]], recursive=TRUE) ], na.rm=TRUE)
#    x.max = max(x1, na.rm=TRUE)
#    x1[ x1 < 0 ] = 0
#    x1[ x1 >= x.max ] = x.max
# print(range(x1, na.rm=TRUE))
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
  if (g %in% rownames(ep.prediction)) {

    # plot prediction
    x1 = ep.prediction[g,]

#    x.max = max(x1[ c(lin.list[[root]], recursive=TRUE) ], na.rm=TRUE)
    x.max = max(x1, na.rm=TRUE)         # using the overall maximum for now
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
    plot.expr.in.fractions(m.cell, x.fraction[g,], reporters, root, read.counts[g,])
  }
}

plot.it = function() {

  pdf(paste(output.dir, "/predictions_P0.pdf", sep=""), width=12, height=12.5)
  par(mfrow=c(4,1))
  par(mar=c(0.1,8,3,3)+0.1)
  plot.prediction.and.fractions.all(abplp.reporters, m.cell.orig, r.m.scaled,
    expr.cell, trunc.pseudo.p, "P0", r1)
  dev.off()
return()
  pdf(paste(output.dir, "/predictions_ABplp.pdf", sep=""), width=5, height=12.5)
  par(mfrow=c(4,1))
  par(mar=c(0.1,8,3,3)+0.1)
  plot.prediction.and.fractions.all(abplp.reporters, m.cell.orig, r.m.scaled,
    expr.cell, trunc.pseudo.p, "ABplp", r1)
  dev.off()

  pdf(paste(output.dir, "/predictions_ABprp.pdf", sep=""), width=5, height=12.5)
  par(mfrow=c(4,1))
  par(mar=c(0.1,8,3,3)+0.1)
  plot.prediction.and.fractions.all(abplp.reporters, m.cell.orig, r.m.scaled,
    expr.cell, trunc.pseudo.p, "ABprp", r1)
  dev.off()
}

