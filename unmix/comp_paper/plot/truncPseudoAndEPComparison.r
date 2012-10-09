# Comparison of the truncated pseudoinverse and EP.

library("clinfun")

source("R/unmix/eval.r")
# source("R/unmix/comp_paper/plot/load.unmixed.results.2.r")

x.tp = unmix.r[["tp"]][["expr.cell"]][[30]]
x.ep = unmix.r[["ep"]][["expr.cell"]][[30]]

a = data.frame(
  tp.cor = diag(cor(t(expr.cell), t(x.tp))),
  ep.cor = diag(cor(t(expr.cell), t(x.ep))),
  tp.auc = auc(expr.on.off[["expr.cell"]], x.tp),
  ep.auc = auc(expr.on.off[["expr.cell"]], x.ep))

plot.it = function() {
  pdf("git/unmix/comp_paper/plot/truncPseudoAndEPComparison.pdf",
    width=10, height=5)
  par(mfrow=c(1,2))

  plot(a$tp.cor, a$ep.cor,
    main="Correlation",
    xlab="Constrained pseudoinverse", ylab="Expectation propagation",
    pch=20, cex=0.5)
  g = rownames(a)[ abs(a$tp.cor - a$ep.cor) >= 0.1 ]
  text(a[g,"tp.cor"], a[g,"ep.cor"], labels=g, cex=0.65, pos=4, offset=0.2)

  plot(a$tp.auc, a$ep.auc,
    main="Area under curve",
    xlab="Constrained pseudoinverse", ylab="Expectation propagation",
    pch=20, cex=0.5)
  g = rownames(a)[ abs(a$tp.auc - a$ep.auc) >= 0.1 ]
  text(a[g,"tp.auc"], a[g,"ep.auc"], labels=g, cex=0.65, pos=4, offset=0.2)

  dev.off()
}

plot.it()

