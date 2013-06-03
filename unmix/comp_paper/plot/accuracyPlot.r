
library("clinfun")

source("R/unmix/eval.r")

# source("R/unmix/comp_paper/plot/load.unmixed.results.2.r")

xlim = c(0, max(num.reporters))

#method.name = c("ep", "tpc", "tp", "pseudo")
#method.label = c("expectation propagation", "pseudoinverse w/ correlation",
#  "pseudoinverse", "naive pseudoinverse")
#method.color = c("black", "red", "green", "blue")

method.name = c("ep", "tp", "tpc", "pseudo", "mf")
method.label = c("expectation propagation",
  "constrained pseudoinverse",
  "constr. pseudoinverse w/ correlation",
  "pseudoinverse", "average fold-change")
method.color = c("black", "red", "purple", "blue", "green")
# method.color = c("#80808080", "#80000080", "#00800080", "#00008080")

# the "actual" answers
expr = list(expr.cell = expr.cell,
  sim.cell.cor = sim.cell.cor,
  sim.expr.cell = sim.expr.cell,
  sim.expr.sym = sim.expr.sym)

# definition of when genes are "on" or "off" (somewhat hokey)
expr.on.off = list(expr.cell = expr.cell >= 2000,
  sim.cell.cor = sim.cell.cor >= 0.2,
  sim.expr.cell = sim.expr.cell >= 5,
  sim.expr.sym = sim.expr.sym >= 5)

# Plots "correlation" accuracy for several data sets w.r.t. number of reporters.
plot.correlation.accuracy = function() {

  # way to compute accuracy
  accuracy = function(x, x.d)
    sapply(num.reporters, function(n) mean(cor1(x, x.d[[n]])))

  plot.it = function(data.set, data.set.name) {
    plot(0,0, xlim=xlim, ylim=c(0,1), type="n",
      main=data.set.name, xlab="Number of reporters", ylab="Correlation",
      cex.main=1.1, cex.axis=1, cex.lab=1)
    if (data.set == "expr.cell") {
      legend("bottomright", legend=method.label, col=method.color,
        lwd=1, seg.len=2)
#      mtext("a)", adj=0, cex=2)
    }

    for(i in 1:5) {
      par(new=TRUE)
      a = accuracy(expr[[data.set]], unmix.r[[method.name[i] ]][[data.set]])
      plot(num.reporters-(i%%2), a, type="b",
        xlim=xlim, ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n",
        col=method.color[[i]], pch=20, cex=0.6)
    }
  }

  plot.it("expr.cell", "Measured expression")
  plot.it("sim.cell.cor", "Simulated patterns based on correlation")
  plot.it("sim.expr.cell", "Simulated one-lineage patterns")
  plot.it("sim.expr.sym", "Simulated two-lineage symmetric patterns")
}

# Plots "area under curve" accuracy for several data sets
# w.r.t. number of reporters.
plot.roc.accuracy = function() {

  # Computes AUC.
  compute.auc = function(x.on.off, x.prediction) {
    x.on.off = ifelse(as.matrix(x.on.off) >= 0.5, 1, 0)

    n = dim(x.on.off)[1]
    auc = rep(NA, n)

    for(i in 1:n) {
      if (sd(x.on.off[i,]) > 0 && sum(is.na(x.prediction[i,])==0)) {
#        cat(i, "")
        a = roc.area.test(x.prediction[i,], x.on.off[i,])
        auc[i] = as.numeric(a$area)
      }
    }

    mean(auc, na.rm=TRUE)
  }

  # Accuracy by AUC.
  accuracy = function(x, x.d)
    sapply(num.reporters, function(n) compute.auc(x, x.d[[n]]))
  ylim = c(0.75, 1)

  plot.it = function(data.set, data.set.name) {
    plot(0,0, xlim=xlim, ylim=ylim, type="n",
      main=data.set.name, xlab="Number of reporters", ylab="Area under curve",
      cex.main=1.1, cex.axis=1, cex.lab=1)
    if (data.set == "expr.cell") {
      legend("bottomright", legend=method.label, col=method.color,
        lwd=1, seg.len=2)
#      mtext("b)", outer=TRUE, cex=2)
    }

    for(i in 1:5) {
      par(new=TRUE)
      a = accuracy(expr.on.off[[data.set]], unmix.r[[method.name[i] ]][[data.set]])
      plot(num.reporters-(i%%2), a, type="b",
        xlim=xlim, ylim=ylim, xlab="", ylab="", xaxt="n", yaxt="n",
        col=method.color[[i]], pch=20, cex=0.6)
    }
  }

  plot.it("expr.cell", "Measured expression")
  plot.it("sim.cell.cor", "Simulated patterns based on correlation")
  plot.it("sim.expr.cell", "Simulated one-lineage patterns")
  plot.it("sim.expr.sym", "Simulated two-lineage symmetric patterns")
}

plot.it = function() {
  pdf("git/unmix/comp_paper/plot/accuracyPlot.pdf",
    title="Accuracy vs. number of reporters", width=11.5, height=7)
  par(mfrow=c(2,4))

  plot.correlation.accuracy()

  plot.roc.accuracy()

  dev.off()
}

plot.it()

