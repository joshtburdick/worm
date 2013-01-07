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

source("multiple_starting_points.r")

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

# Summarizes samples at several windows.
summarize.samples.1 = function(x, window.size) {
  n = floor( nrow(x) / window.size )
  r = array(dim=c(n, 5, ncol(x)),
    dimnames=list(window=c(1:n), stat=NULL, cell=NULL))

  for(i in 1:n) {
    i1 = (i-1) * window.size + 1
    i2 = i * window.size
cat(i1,i2)
    s = summarize.samples(x[c(i1:i2),])
    r[i,,] = s
    dimnames(r)[[2]] = rownames(s)
    dimnames(r)[[3]] = colnames(s)
  }

  r
}

# now includes several restarts
unmix.xsample = function(m, x.fraction) {
  x.summary = NULL
  X1 = NULL
  for(iter in 1:3) {
    X = run.xsample(m, x.fraction, 5000, 50000)
    x.summary[[iter]] = summarize.samples.1(X, 5000)
    X1[[iter]] = X
#    x.summary[[iter]] <- summarize.samples(X)
  }

  list(x = as.vector(apply(x.summary[[1]][,"mean",], 2, mean)),
    x.summary = x.summary, X = X1)
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

unmix.xsample.1 = function(gene, f) {
cat("file = ", f, "\n")
  unmix.result = NULL
  st = system.time( unmix.result <-
    run.unmix.1(expr.cell[gene,,drop=FALSE], m.cell, unmix.xsample, reporters$picked, 30) )
  unmix.result$system.time = st

  save(unmix.result, file=f)
}

# Runs unmix.xsample.1, several times.
unmix.xsample.2 = function(gene) {
  outdir = "multiple_restart_pseudoinverse_start"
  system(paste("mkdir -p ", outdir))

  unmix.xsample.1(gene, paste(outdir, "/", gene, ".sampling.Rdata",
    sep="", collapse=""))
}

a = commandArgs(trailingOnly=TRUE)[-1]
gene = sub(".sampling", "", a[1])

unmix.xsample.2(gene)

