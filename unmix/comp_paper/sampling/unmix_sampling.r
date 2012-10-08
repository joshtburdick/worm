# Does unmixing by sampling.
# XXX bit of a hack

library("corpcor")
library("limSolve")

source("~/gcb/R/unmix/eval.r")
source("R/unmix/comp_paper/sampling/xsample1.r")

load("~/gcb/git/unmix/unmix_comp/data/tree_utils.Rdata")

wd = getwd()
setwd("~/gcb/git/unmix/unmix_comp/src/")
source("unmix_test.r")
source("pseudoinverse/unmix.r")
setwd(wd)

num.cells = apply(cell.lineage.matrix, 1, sum)


# Does unmixing of one gene using sampling.
# Args:
#   m - sort matrix
#   x.fraction - amounts of expression in each fraction
# Returns: function which, given a sort matrix and
#   expression in each fraction, returns unmixed expression.
unmix.xsample = function(m, x.fraction) {
  num.cells = dim(m)[2]

  r = xsample1(E=m, F=as.vector(x.fraction),
    G=diag(num.cells), H=rep(0, num.cells),
    tol=1e-3, type="rda",
    burninlength=5000, iter=50000, test=FALSE)

  x = as.vector(apply(r$X, 2, mean))
  names(x) = colnames(m)
  list(x = x, sampling.x = r$X)
}

unmix.xsample.1 = function() {
  nr = 30

  unmix.result = NULL
  st = system.time( unmix.result <-
    run.unmix.1(expr.cell, m.cell, unmix.xsample, reporters$picked, nr) )
  unmix.result$system.time = st
  unmix.result
  save(unmix.result, file="git/unmix/comp_paper/unmix_sampling_1.Rdata")
}


