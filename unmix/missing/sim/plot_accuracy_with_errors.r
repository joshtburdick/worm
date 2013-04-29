# Plots accuracy, when several types of errors are present.

source("git/plot/label_panel.r")
source("git/unmix/missing/sim/plot_accuracy_with_noise.r")

load("git/unmix/missing/sim/accuracy.missing.cells.Rdata")
load("git/unmix/missing/sim/accuracy.add.Rdata")
load("git/unmix/missing/sim/accuracy.delete.Rdata")
load("git/unmix/missing/sim/accuracy.flip.Rdata")

# Adds a column to a result dataset indicating
# which run it came from (arguably, this should have
# been in the data set in the first place.)
group.results = function(r) {
  n = nrow(r) / 123
  stopifnot(n == trunc(n))
  r = cbind(r, g = rep(1:n, each=123))
  r
}

# Summarizes accuracy of a set of results.
summarize.results = function(r) {
  r = group.results(r)

  a = data.frame(
    num.cells.diff =
      c(by(r$num.cells.diff, r$g, function(x) mean(x, na.rm=TRUE))),
    cor.incorrect.m =
      c(by(r$cor.incorrect.m, r$g, function(x) mean(x, na.rm=TRUE))),
    auc.incorrect.m =
      c(by(r$auc.incorrect.m, r$g, function(x) mean(x, na.rm=TRUE))))
  a
}

summary.missing.cells = summarize.results(accuracy.missing.cells)
summary.add = summarize.results(accuracy.add)
summary.delete = summarize.results(accuracy.delete)
summary.flip = summarize.results(accuracy.flip)

# Plots accuracy with various numbers of errors.
plot.accuracy.with.errors = function() {
  pdf("git/unmix/missing/sim/accuracyWithErrors.pdf",
    width=6, height=7)
  par(mfrow=c(3,2))
  par(mar=c(4,4,4,1) + 0.1)

  entries.in.matrix = 30 * 1341

  plot(100*summary.flip$num.cells.diff / entries.in.matrix,
    summary.flip$auc.incorrect.m, xlim=c(0,5), # c(0,2000),
    ylim=c(0.6, 1), pch=20, cex=0.75, col="#00000080",
    main="Accuracy with sort errors",
    xlab="Percentage of entries changed", ylab="Area under curve")
  label.panel("a)", gp=gpar(fontsize=16))
  plot(100*summary.flip$num.cells.diff / entries.in.matrix,
    summary.flip$cor.incorrect.m, xlim=c(0,5),
    ylim=c(0, 1), pch=20, cex=0.75, col="#00000080",
    main="Accuracy with sort errors",
    xlab="Percentage of entries changed", ylab="Correlation")

  plot(100*summary.missing.cells$num.cells.diff / entries.in.matrix,
    summary.missing.cells$auc.incorrect.m, xlim=c(0,5),
    ylim=c(0.6, 1), pch=20, cex=0.75, col="#00000080",
    main="Accuracy with missing cells",
    xlab="Percentage of entries changed", ylab="Area under curve")
  label.panel("b)", gp=gpar(fontsize=16))
  plot(100*summary.missing.cells$num.cells.diff / entries.in.matrix,
    summary.missing.cells$cor.incorrect.m, xlim=c(0,5),
    ylim=c(0, 1), pch=20, cex=0.75, col="#00000080",
    main="Accuracy with missing cells",
    xlab="Percentage of entries changed", ylab="Correlation")

if (FALSE) {
  plot(summary.delete$num.cells.diff,
    summary.delete$auc.incorrect.m,
    ylim=c(0.6, 1), pch=20, cex=0.75, col="#00000080",
    main="Accuracy with false positives",
    xlab="Number of entries changed", ylab="Area under curve")
  label.panel("c)", gp=gpar(fontsize=16))
  plot(summary.delete$num.cells.diff,
    summary.delete$cor.incorrect.m,
    ylim=c(0, 1), pch=20, cex=0.75, col="#00000080",
    main="Accuracy with false positives",
    xlab="Number of entries changed", ylab="Correlation")

  plot(summary.add$num.cells.diff,
    summary.add$auc.incorrect.m,
    ylim=c(0.6, 1), pch=20, cex=0.75, col="#00000080",
    main="Accuracy with false negatives",
    xlab="Number of entries changed", ylab="Area under curve")
  label.panel("d)", gp=gpar(fontsize=16))
  plot(summary.add$num.cells.diff,
    summary.add$cor.incorrect.m,
    ylim=c(0, 1), pch=20, cex=0.75, col="#00000080",
    main="Accuracy with false negatives",
    xlab="Number of entries changed", ylab="Correlation")
}

  plot.accuracy.with.measurement.noise()

  dev.off()
}

# Does regression to measure error with different numbers
# of incorrect cells.
incorrect.m.regression = function() {
  r = summary.flip
  r$fraction.diff = r$num.cells.diff / (30 * 1341)
  r = r[ r$fraction.diff < 0.05 , ]

  r$log.auc = log(r$auc.incorrect.m)
  m.auc = lm(log.auc ~ fraction.diff, data=r)
  predicted.auc = exp(predict.lm(m.auc, newdata=data.frame(fraction.diff = c(0, 0.05))))
  cat("AUC:\n")
  print(c(predicted.auc, (predicted.auc[1] - predicted.auc[2]) / predicted.auc[1]))

  r$log.cor = log(r$cor.incorrect.m)
  m.cor = lm(log.cor ~ fraction.diff, data=r)
  predicted.cor = exp(predict.lm(m.cor, newdata=data.frame(fraction.diff = c(0, 0.05))))
  cat("correlation:\n")
  print(c(predicted.cor, (predicted.cor[1] - predicted.cor[2]) / predicted.cor[1]))
}

plot.accuracy.with.errors()

