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
    width=10, height=10)
  par(mfrow=c(2,2))

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

  plot(a$tp.cor, a$tp.auc,
    main="Constrained pseudoinverse",
    xlab="Correlation", ylab="Area under curve",
    pch=20, cex=0.5)
  g = unique(c(rownames(a)[a$tp.cor < 0.7 | a$tp.auc < 0.85],
    sample(rownames(a), 8)))
  text(a[g,"tp.cor"], a[g,"tp.auc"], labels=g, cex=0.65, pos=4, offset=0.2)

  plot(a$ep.cor, a$ep.auc,
    main="Expectation propagation",
    xlab="Correlation", ylab="Area under curve",
    pch=20, cex=0.5)
  g = unique(c(rownames(a)[a$ep.cor < 0.74 | a$ep.auc < 0.88],
    sample(rownames(a), 8)))
  text(a[g,"ep.cor"], a[g,"ep.auc"], labels=g, cex=0.65, pos=4, offset=0.2)
  dev.off()
}

plot.it()

