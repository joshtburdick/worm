# Does unmixing by sampling.
# XXX bit of a hack

library("corpcor")
library("limSolve")

wd = getwd()

setwd("~/gcb/")
source("R/unmix/eval.r")

# source("R/unmix/comp_paper/sampling/xsample1.r")
setwd("~/gcb/git/unmix/unmix_comp/src/")
source("sampling/cdaCpp.r")

setwd(wd)

load("~/gcb/git/unmix/unmix_comp/data/tree_utils.Rdata")

wd = getwd()
setwd("~/gcb/git/unmix/unmix_comp/src/")
source("unmix_test.r")
source("pseudoinverse/unmix.r")
setwd(wd)

source("multiple_starting_points.r")

num.cells = apply(cell.lineage.matrix, 1, sum)

# Runs xsample() with sensible-seeming parameters.
# This only includes the cells for which lsei() indicates the 
# "most likely estimate" isn't zero.
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
    x0 = overdispersed.start.1(m[,cl], x.fraction, 0),  # 0.9),
#    jmp=1 * max(x.fraction) / (sum(cl)),   # was 10 * that
    burninlength=burninlength, iter=iter, test=FALSE)

  # add zero cells back in
  x = matrix(0, nrow=iter, ncol=num.cells)
  colnames(x) = colnames(m)
  x[,cl] = r$X
  x
}

# Summarizes samples.
# Args:
#   x - matrix of samples
# Returns: matrix of summary statistics
summarize.samples = function(x) {
  r = rbind(apply(x, 2, mean),
    apply(x, 2, var),
    apply(x, 2, median),
    apply(x, 2, min),
    apply(x, 2, max))
  rownames(r) = c("mean", "var", "median", "min", "max")
  r
}

# Does unmixing. Keeps track of statistics in an efficient way,
# allowing doing very long runs with less memory.
# Args:
#   iters - number of samples
#   window.size - number of samples to group together
#     (to reduce size of output statistics)
#   burnin.windows - number of initial windows to discard
#     (thus, number of burn-in samples is this * window.size)
#   m - the matrix defining reporter fractions
#   x.fraction - amount in each reporter
# Returns: list with elements
#   x - the estimate (average of all samples after burn-in)
#   x.summary - array containing summary statistics
#   x.last - the last sample obtained
run.xsample.lowmem = function(iters, window.size, burnin.windows=0)
    function(m, x.fraction) {
  overdisperse.alpha = 0.9

  num.cells = dim(m)[2]

  # XXX hard-coded for now
  thinning = 2000

  # find initial estimate
  x.fraction.orig = x.fraction
  x.fraction = x.fraction
cat("about to call lsei\n")
  x0 = lsei(A=diag(num.cells), B=rep(0, num.cells),
    E=m, F=x.fraction, G=diag(num.cells), H=rep(0, num.cells),
    tol=1e-4, type=2)$X
cat("after lsei\n")

  # find cells estimated to be zero
  cl = which(abs(x0) >= 1e-30)
  cat(length(cl), "were non-zero\n")

  # x0 is the non-zero entries of that estimate
  x0 = x0[cl]

  # possibly use an "overdispersed" starting point
#  if (!is.null(overdisperse.alpha)) {
#    x0 = overdispersed.start.1(m[,cl], x.fraction, overdisperse.alpha)
#  }

  # allocate array for summary stats of samples
  num.windows = trunc( iters / window.size )
  stopifnot(burnin.windows < num.windows)
  iters.per.window = iters / num.windows
  x.summary = array(0, dim=c(num.windows, 5, num.cells),
    dimnames=list(window=NULL,
      stat=c("mean", "var", "median", "min", "max"),
      cell=NULL))

  # do sampling
  for(i in 1:num.windows) {
cat("window =", i, "\n")
    X = sample.cda(m[,cl], x.fraction, x0,
      iters.per.window, thinning)
# print(dim(X))
# print(X[1:4,1:4])
    x.summary[i,,cl] = summarize.samples(X)

    x0 = as.vector(X[ nrow(X) , ])
  }

  # compute average (possibly discarding burn-in)
  x = if (burnin.windows == 0)
    as.vector(apply(x.summary[,"mean",], 2, mean))
  else
    as.vector(apply(x.summary[-c(1:burnin.windows),"mean",], 2, mean))
  names(x) = colnames(m)
  dimnames(x.summary)[[3]] = colnames(m)

  list(x = x, burnin.windows = burnin.windows, x.summary = x.summary )
}

# Does sampling, including several restarts.
unmix.multiple.restarts =
  function(iters, window.size, burnin.windows=0, num.restarts=10)
    function(m, x.fraction) {
  r = list()

  for(iter in 1:num.restarts) {
    r[[iter]] = run.xsample.lowmem(iters, window.size, burnin.windows=0)(m, x.fraction)
  }

  list(x = r[[1]][["x"]], r = r)
}

# Does unmixing for one gene
unmix.xsample.1 = function(gene, f) {
# cat("file = ", f, "\n")
  unmix.result = NULL
  st = system.time( unmix.result <-
    run.unmix.1(expr.cell[gene,,drop=FALSE], m.cell,
#      run.xsample.lowmem(10, 5, 0),
      unmix.multiple.restarts(50000, 1000, 0, 5),
      reporters$picked, 30) )
  unmix.result$system.time = st

  save(unmix.result, file=f)
}

# Runs unmix.xsample.1, several times.
unmix.xsample.2 = function(gene) {
  outdir = "multiple_restart"
  system(paste("mkdir -p ", outdir))

  unmix.xsample.1(gene, paste(outdir, "/", gene, ".Rdata",
    sep="", collapse=""))
}

a = commandArgs(trailingOnly=TRUE)[-1]
gene = sub(".sampling", "", a[1])

unmix.xsample.2(gene)

