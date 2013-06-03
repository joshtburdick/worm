# Plots sampling results.

source("R/plot/plot_expr.r")

load("R/unmix/comp_paper/expr.cell.Rdata")

load("git/unmix/comp_paper/sampling/plot/samplingStats.Rdata")

# Plots a graph of one gene's prediction, and actual expression.
plot.comparison.linear = function(expr.cell, sampling.stats, gene.name) {
  samples = sampling.stats[gene.name,,]

  if(any(is.na(samples))) {
    return()
  }

  y.lo = pmax(0, samples["mean",] - 2 * samples["sd",])
  y.hi = samples["mean",] + 2 * samples["sd",]
  ymax = max(expr.cell[gene.name,], y.lo, y.hi)

  plot(1,0, xlim=c(1,1341), ylim=c(0,ymax), type="n", main=gene.name,
    xlab="", ylab="expression", xaxt="n", cex.axis=0.5, las=1)
  segments(1:1341, y.lo, 1:1341, y.hi, col="lightgrey")
  par(new=TRUE)
  plot(1:1341, as.vector(expr.cell[gene.name,]), ylim=c(0,ymax),
    xaxt="n", yaxt="n", xlab="", ylab="", pch=20, cex=0.2)

#  axis(1, at=cell.to.column[int.n], labels=int.n, las=2, cex.axis=0.35)
  axis(1, at=cell.to.column[lin.12.cell], labels=lin.12.cell,
    cex.axis=0.7)
}

# Plots all the genes.
plot.all = function(sampling.stats, output.dir) {
  options(digits=18)

  system(paste("mkdir -p ", output.dir))

  for(gene.name in dimnames(sampling.stats)[[1]]) {
cat(gene.name, "")
    pdf(paste(output.dir, "/", gene.name, ".pdf", sep="", collapse=""),
      width=9, height=3)
    plot.comparison.linear(expr.cell, sampling.stats, gene.name)
    dev.off()
  }

  options(digits=7)
}

if (TRUE) {
  plot.all(samplingStats,
    "git/unmix/comp_paper/sampling/plot/linear/")
}


