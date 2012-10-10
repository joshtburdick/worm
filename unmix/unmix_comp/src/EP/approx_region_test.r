# Tests of approx.region .

library(limSolve)

source("EP/approx_region.r")

# Compares xsample() with approx.region on some simple cases.
# Args:
#   m, n - number of rows and columns in A
# Returns: list with components:
#   A, x, b - the random matrix, x, and b
#   ep.r - the result from running approx.region()
#   xsample.result - the result from running xsample()
ep.xsample.comparison = function(m, n) {

  A = matrix(rbeta(m*n, 0.5, 0.5), nrow=m, ncol=n)
  x = rbeta(n, 0.5, 0.5)
  b = as.vector( A %*% x )

  xsample.r = xsample(E = A, F = b, G = diag(n), H = rep(0, n),
    iter=1000, burninlength=1000)
  xsample.r$m = apply(xsample.r$X, 2, mean)
  xsample.r$v = apply(xsample.r$X, 2, var)
  ep.r = approx.region(A, b, 0*b)

  list(A = A, x = x, b = b, xsample.r = xsample.r, ep.r = ep.r)
}

# Plots some mean and variance comparisons.
plot.mv.comparisons = function() {
  system.time(r <- ep.xsample.comparison(20,100)) 

  plot(r$xsample.r$m, r$ep.r$m)
  plot(sqrt(r$xsample.r$v), sqrt(r$ep.r$v))
}


