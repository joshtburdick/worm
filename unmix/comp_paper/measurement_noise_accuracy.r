# Measures accuracy with measurement noise added.

library("corpcor")

source("~/gcb/R/unmix/eval.r")

load("~/gcb/git/unmix/unmix_comp/data/tree_utils.Rdata")

wd = getwd()
setwd("~/gcb/git/unmix/unmix_comp/src/")
source("unmix_test.r")
# source("pseudoinverse/unmix.r")
source("EP/unmix.r")
setwd(wd)

num.cells = apply(cell.lineage.matrix, 1, sum)

# Runs an unmixing function with different numbers of reporters
# (cross-validated), and with noise added to the measurements.
# Args:
#   x - the expression matrix
#   m - the sorting matrix
#   unmix.f - the unmixing function
#   reporter.list - the reporters
#   nr - number of reporters
#   noise.sd - the amount of noise: noise will be obtained by
#     multiplying the fraction amounts by Normal(1, noise.sd)
run.unmix.perturbed.matrix = function(x, m, unmix.f, reporter.list, nr, noise.sd) {
  x.predicted = matrix(nrow=nrow(x), ncol=ncol(x))
  dimnames(x.predicted) = dimnames(x)

  for(g in rownames(x)) {
 cat(g, "")

    # set of reporters, with another one added in if g is present
    r1 = if (g %in% reporter.list[1:nr])
      setdiff(reporter.list[1:(nr+1)], g)
    else
      reporter.list[1:nr]
# cat("r1 =", r1, "\n")

    m1 = m[ c("all", r1), ]
    x.fraction = x[g,] %*% t(m1)
    x.fraction = x.fraction *
      rnorm(nrow(x.fraction), mean = 1, sd = noise.sd)
    try({
      r = unmix.f(m1, as.vector(x.fraction))
      x.predicted[g,] = as.vector(r$x)
    })
  }

  rownames(x.predicted) = rownames(x)

  x.predicted
}

# Computes accuracy for some noise level.
compute.accuracy.with.noise = function(noise.sd) {
cat(noise.sd, "\n")
  x.unmix = run.unmix.perturbed.matrix(expr.cell, m.cell, unmix.ep,
    reporters$picked, 30, noise.sd)

  c(noise.sd = noise.sd,
    cor = mean(diag(cor(t(expr.cell), t(x.unmix))), na.rm=TRUE),
    auc = mean(auc(expr.cell >= 2000, x.unmix), na.rm=TRUE),
    num.non.na = sum(apply(!is.na(x.unmix), 1, sum) > 0))
}

accuracy.with.noise = NULL
for(iter in 1:3)
  accuracy.with.noise = rbind(accuracy.with.noise,
    t( sapply(c(0:20)/20, compute.accuracy.with.noise)))

write.table(accuracy.with.noise,
  file="git/unmix/comp_paper/measurement_noise_accuracy.tsv",
  col.names=NA, sep="\t")

plot.it = function() {
  pdf("git/unmix/comp_paper/accuracy with measurement noise.pdf",
    width=9, height=4)
  par(mfrow=c(1,2))

  plot(accuracy.with.noise[,"noise.sd"], accuracy.with.noise[,"cor"],
    ylim=c(0.5,1), pch=20, cex=0.75,
    xlab="noise standard deviation", ylab="correlation",
    main="Accuracy with measurement noise")

  plot(accuracy.with.noise[,"noise.sd"], accuracy.with.noise[,"auc"],
    ylim=c(0.5,1), pch=20, cex=0.75,
    xlab="noise standard deviation", ylab="area under the curve",
    main="Accuracy with measurement noise")

  dev.off()
}

#  write.table(missing.cells.separate.fractions, col.names=NA,
#    file="git/unmix/missing/sim/missing_cells_separate_fractions.txt", sep="\t")

# foo = run.unmix.perturbed.matrix(expr.cell, m.cell, unmix.pseudoinverse,
#   reporters$picked, 30, 1)


plot.it()

