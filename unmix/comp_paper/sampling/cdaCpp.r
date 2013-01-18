# Attempt at implementing somthing like the "cda"
# variant of xsample()'s hit-and-run sampling method.

library("Rcpp")
library("limSolve")

sourceCpp(file = "git/unmix/comp_paper/sampling/cda.cpp")

# Does random directions sampling, jumping in one of
# a set of directions (which presumably span the nullspace.)
# Args:
#   A, b - these define the constraint that Ax = b
#   x0 - the starting point
#   num.samples - the number of samples to include
# not yet implemented:
#   thinning - number of samples to skip (e.g., setting this
#     to 1000 will only store every thousandth sample.)
#     Useful because this method has high autocorrelation.
#     By default, set to 1 (which returns every sample.)
#   renormalization interval?
# Returns: matrix of samples
sample.cda = function(A, b, x0, num.samples, thinning) {

  # FIXME compute x0?

  # compute the nullspace
  Z = t(Null(t(A)))

  r = cdaCppCore(A, b, x0, Z, num.samples, thinning)
  r
}

# Tiny test; samples three numbers adding up to one.
cdaCpp.test0 = function() {
  A = matrix(c(1,1,1)/3, ncol=3, byrow=TRUE)
  b = c(1)
#  x0 = c(1,1,1) / 3
  x0 = c(0.99,0.005,0.005)
#  Z = matrix(c(2,-1,-1, -1,2,-1, -1,-1,2), ncol=3, byrow=TRUE)
  Z = matrix(c(2,-1,-1, -1,2,-1), ncol=3, byrow=TRUE)
  r = cdaCppCore(A, b, x0, Z, 100000, 10)

  png("git/unmix/comp_paper/sampling/cdaCppTest0.png",
    width=900, height=300)
  par(mfrow=c(1,3))
  for(j in c(1:3))
    hist(r[,j], breaks=100, col="grey", xlim=c(0,1))
  dev.off()
}

