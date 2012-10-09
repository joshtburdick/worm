# Plots some actual expression patterns, and predictions
# in a sub-lineage.

# source("R/unmix/comp_paper/plot/load.unmixed.results.2.r")
# source("git/plot/plot_expr.r")

expr.prediction = unmix.r[["pseudo"]][["expr.cell"]][[30]][,lin.node.names]

# Plots part of the C lineage.
# Args:
#   expr - the expression pattern
plot.expr.1 = function(expr, main, max.intensity, max.color) {
  par(mar=c(2,4,1.8,2)+0.1)
  expr[expr > max.intensity] = max.intensity

  c1 = as.vector(col2rgb(max.color) / 255)
  r = rgb(scale.to.unit(expr) * c1[1],
    scale.to.unit(expr) * c1[2],
    scale.to.unit(expr) * c1[3])
  names(r) = names(expr)
  plot.segments.per.cell(r, main, root="C", lwd=3)
  title(main)
}

fig.1 = function() {
  pdf("git/unmix/comp_paper/plot/1. illustration of method.pdf",
    width=8, height=8)
  par(mfcol=c(3,2))

  plot.expr.1(expr.cell["alr-1",], "alr-1", 2000, "red")
  plot.expr.1(expr.cell["nhr-171",], "nhr-171", 2000, "orange")
  plot.expr.1(expr.cell["tlp-1",], "tlp-1", 3000, "green")

  plot.expr.1(expr.cell["pes-1",], "pes-1", 1400, "purple")
  plot.new()
  plot.expr.1(expr.prediction["pes-1",], "pes-1 prediction", 1450, "purple")

  dev.off()
}

fig.1()


