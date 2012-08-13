
# Plots cross-validated unmixing predictions.
plot.it = function() {
  output.dir = "git/unmix/eval/comparison/"
  system(paste("mkdir -p", output.dir))

  for(g in reporters) {
    png(paste(output.dir, "/", g, ".png", sep=""), width=1600, height=1200)
    par(mfrow=c(4,1))
    par(mar=c(2,4,4,2)+0.1)

    a = rgb(scale.to.unit(expr.cell[g,]), 0, 0)
    names(a) = colnames(expr.cell)
    plot.segments.per.cell(a, paste(g, " from imaging", sep=""),
      root="P0", times=c(-20,550), lwd=3)

    a = rgb(0, scale.to.unit(r.pseudo$x[g,] / M["all",]), 0)
    names(a) = colnames(M)
    plot.segments.per.cell(a, paste(g, " pseudoinverse (not using ", g, ")", sep=""),
      root="P0", times=c(-20,550), lwd=3)

    a = rgb(0, scale.to.unit(r.lsei$x[g,] / M["all",]), 0)
    names(a) = colnames(M)
    plot.segments.per.cell(a, paste(g, " pseudoinverse > 0 (not using ", g, ")", sep=""),
      root="P0", times=c(-20,550), lwd=3)

    a = rgb(0, scale.to.unit(r.ep$x[g,] * M["all",]), 0)
    names(a) = colnames(M)
    plot.segments.per.cell(a, paste(g, " EP (not using ", g, ")", sep=""),
      root="P0", times=c(-20,550), lwd=3)

    dev.off()
  }
}


