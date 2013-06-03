# Computes the mean and variance for each sublineage of cells.

library("Matrix")

source("unmix_test.r")   # for loading the cell-sorting matrix
source("EP/approx_region.r")

load("../data/tree_utils.Rdata")

num.cells = apply(cell.lineage.matrix, 1, sum)

# Given an EP result, gets the mean and variance for each sublineage.
# Args:
#   m - the cell-sorting matrix
#   ep - an EP result, including "t" and "reporters"
#   original.scale - should this include the original scaling?
# Returns: list with elements:
#   m - the mean
#   V - the full covariance matrix
#   lineage.m, lineage.v - mean and variance of lineage totals
get.lineage.totals.1 = function(m, ep, original.scale) {

  f = function(a) {
    a[ is.na(a) ] = 0
    a[ a < 0 ] = 0
    a
  }

  # get mean and covariance, by adding in the linear constraint
  m1 = m[ep$reporters,]
  r = lin.constraint(ep$t[,"m"], ep$t[,"v"], m1, ep$x.f, 0 * ep$x.f)

  s = if (original.scale) ep$scale else 1

  # restore the original scale (and subtract off epsilon)
  per.cell.m = s * (r$m - ep$eps)

  per.cell.V = s * s * r$V

  a = list(per.cell.mean = f(per.cell.m),
    per.cell.var = f(diag(per.cell.V)),
    lineage.mean=
      f(cell.lineage.matrix %*% per.cell.m / num.cells),
    lineage.var =
      f(diag(cell.lineage.matrix %*% per.cell.V %*% t(cell.lineage.matrix)) /
      num.cells))
}

# Gets mean and variance for each gene.
get.lineage.totals = function(file, original.scale) {
  load(file)

  a = array(dim=c(length(unmix.result$r), length(lin.node.names), 4),
    dimnames=c("movie", "cell", "stat"))
  dimnames(a)[[1]] = names(unmix.result$r)
  dimnames(a)[[2]] = lin.node.names
  dimnames(a)[[3]] =
    c("per.cell.mean", "per.cell.var", "lineage.mean", "lineage.var")

  for(g in names(unmix.result$r)) {
    cat(g, "")
    x = get.lineage.totals.1(m.cell, unmix.result$r[[g]], original.scale)
    a[g,,"per.cell.mean"] = as.vector(x$per.cell.mean)
    a[g,,"per.cell.var"] = as.vector(x$per.cell.var)
    a[g,,"lineage.mean"] = as.vector(x$lineage.mean)
    a[g,,"lineage.var"] = as.vector(x$lineage.var)
  }

  a
}

# get lineage totals for one particular case
if (FALSE) {
  lineage.totals = get.lineage.totals("EP/expr.cell.30.Rdata", TRUE)
  save(lineage.totals, file="EP/lineage.totals.Rdata")
  lineage.totals.unscaled = get.lineage.totals("EP/expr.cell.30.Rdata", FALSE)
  save(lineage.totals.unscaled, file="EP/lineage.totals.unscaled.Rdata")
}

lineage.totals = get.lineage.totals("EP.2/expr.cell.30.Rdata", TRUE)
save(lineage.totals, file="EP.2/lineage.totals.Rdata")
lineage.totals.unscaled = get.lineage.totals("EP.2/expr.cell.30.Rdata", FALSE)
save(lineage.totals.unscaled, file="EP.2/lineage.totals.unscaled.Rdata")


