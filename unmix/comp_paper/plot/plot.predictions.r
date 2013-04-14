
source("R/unmix/comp_paper/plot/load.unmixed.results.2.r")
# source("R/plot/plot_expr.r")
source("git/plot/plot_expr.r")
source("git/plot/label_panel.r")


# Plots "actual" expression in red, and two different predictions in green;
# used for comparing methods.
# Args:
#   expr.ref, expr.1 - expression datasets, with one column per cell.
#   output.dir - directory in which to write images
plot.expr.comparison.rgb = function(expr.cell, d.expr.cell, output.dir) {
  system(paste("mkdir -p ", output.dir, sep=""))

  for(g in rownames(expr.cell)) {
cat(g, "")

    pdf(file=paste(output.dir, "/", g, ".pdf", sep=""),
      width=11, height=10)
    par(mar=c(2,4,1.8,2)+0.1)
    par(mfrow=c(3,1))

    # actual expression
    col = rgb(scale.to.unit(expr.cell[g,]),
      scale.to.unit(0*expr.cell[g,]),
      scale.to.unit(0*expr.cell[g,]))
    names(col) = names(expr.cell[g,])
    plot.segments.per.cell(col,
      paste(g, "expression"), lwd=1.3)
    label.panel("a)")


    # prediction with 20 reporters
    col = rgb(scale.to.unit(0*expr.cell[g,]),
      scale.to.unit(d.expr.cell[[20]][g,]),
      scale.to.unit(0*expr.cell[g,]))
    names(col) = names(expr.cell[g,])
    plot.segments.per.cell(col,
      paste(g, "unmixing prediction with 20 reporters"), lwd=1.3)
    label.panel("b)")

    # ... and with 30 reporters
    col = rgb(scale.to.unit(0*expr.cell[g,]),
      scale.to.unit(d.expr.cell[[30]][g,]),
      scale.to.unit(0*expr.cell[g,]))
    names(col) = names(expr.cell[g,])
    plot.segments.per.cell(col,
      paste(g, "unmixing prediction with 30 reporters"), lwd=1.3)
    label.panel("c)")

    dev.off()
  }

cat("\n")
}


plot.expr.comparison.rgb(expr.cell, unmix.r[["tp"]][["expr.cell"]],
  "git/unmix/comp_paper/plot/unmix.comparison/")



