# Does unmixing using lsei().

library("corpcor")

source("unmix_test.r")

# calculate covariance of movies
expr.cell = as.matrix(read.table("../data/exprCell.tsv",
  sep="\t", header=TRUE, row.names=1, as.is=TRUE))
cor.cell.shrunk = cor.shrink(expr.cell, lambda=0.05)

# Does unmixing of one gene.
# Args:
#   si - the covariance matrix
# Returns: function which, given a sort matrix and
#   expression in each fraction, returns unmixed expression.
unmix.lsei = function(si) {
  A = chol(chol2inv(si))

  function(m, x.fraction) {
    num.cells = dim(m)[2]

    # constraints to force everything to be positive
    G = diag(num.cells)
    H = rep(0, num.cells)

    r = lsei(A=A, B=rep(0, num.cells), E=m, F=x.fraction, G=G, H=H, type=1)
    x = as.vector(r$X)
    names(x) = colnames(m)
    list(x = x)
  }
}

run.unmix(unmix.lsei(cor.cell.shrunk), "trunc.pseudoinverse.cor/")

