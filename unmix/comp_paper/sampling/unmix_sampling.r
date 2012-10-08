# Does unmixing by sampling.
# XXX bit of a hack

library("corpcor")
library("limSolve")

wd = getwd()
setwd("~/gcb/")
source("~/gcb/R/unmix/eval.r")
source("R/unmix/comp_paper/sampling/xsample1.r")

load("~/gcb/git/unmix/unmix_comp/data/tree_utils.Rdata")
setwd(wd)

wd = getwd()
setwd("~/gcb/git/unmix/unmix_comp/src/")
source("unmix_test.r")
source("pseudoinverse/unmix.r")
setwd(wd)

num.cells = apply(cell.lineage.matrix, 1, sum)



# Runs xsample() with sensible-seeming parameters.
# This only includes the cells for which lsei() indicates the 
# "most likely estimate" isn't zero
# Args:
#   m - the matrix defining reporter fractions
#   x.fraction - amount in each reporter
#   burninlength, iter - amount of burn-in, and number of iterations
run.xsample = function(m, x.fraction, burninlength, iter) {
  num.cells = dim(m)[2]

  x.fraction.orig = x.fraction
  x.fraction = x.fraction
cat("about to call lsei\n")
  x0 = lsei(A=diag(num.cells), B=rep(0, num.cells),
    E=m, F=x.fraction, G=diag(num.cells), H=rep(0, num.cells),
    tol=1e-4)$X
cat("after lsei\n")

  # find cells estimated to be zero
  cl = which(abs(x0) >= 1e-30)
  cat(length(cl), "were non-zero\n")

  # estimate other cells
  r = xsample1(E=m[,cl], F=x.fraction, G=diag(length(cl)), H=rep(0, length(cl)),
    tol=1e-4, type="rda",
    burninlength=burninlength, iter=iter, test=FALSE)

  # add zero cells back in
  x = matrix(0, nrow=iter, ncol=num.cells)
  colnames(x) = colnames(m)
  x[,cl] = r$X
  x
}

unmix.xsample = function(m, x.fraction) {
  X = run.xsample(m, x.fraction, 5000, 50000)
  r = list(x = apply(X, 2, mean), x.var = apply(X, 2, var),
    x.sampling = X)
  r
}

# Does unmixing of one gene using sampling.
# Args:
#   m - sort matrix
#   x.fraction - amounts of expression in each fraction
# Returns: function which, given a sort matrix and
#   expression in each fraction, returns unmixed expression.
unmix.xsample.old = function(m, x.fraction) {
  num.cells = dim(m)[2]

  r = xsample1(E=m, F=as.vector(x.fraction),
    G=diag(num.cells), H=rep(0, num.cells),
    tol=1e-3, type="rda",
    burninlength=100, iter=100, test=FALSE)

  x = as.vector(apply(r$X, 2, mean))
  names(x) = colnames(m)
  list(x = x, sampling.x = r$X)
}

unmix.xsample.1 = function(gene) {
  unmix.result = NULL
  st = system.time( unmix.result <-
    run.unmix.1(expr.cell[gene,,drop=FALSE], m.cell, unmix.xsample, reporters$picked, 30) )
  unmix.result$system.time = st

  save(unmix.result, file=paste(gene, ".sampling.Rdata", sep=""))
}

a = commandArgs(trailingOnly=TRUE)[-1]
gene = sub(".sampling.Rdata", "", a[1])

unmix.xsample.1(gene)

