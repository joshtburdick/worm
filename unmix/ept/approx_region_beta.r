# Attempts to approximate a region using beta distributions
# (instead of Gaussians.) It's a fairly conservative design,
# in that other than that, it's basically the same.

source("git/unmix/ept/approx_region.r")

# Natural parameters of a beta distribution.
beta.params.to.natural = function(p) {
  x = p - 1
  dimnames(x)[[1]] = c("e1", "e2")
  x
}

# ... and natural params to beta params.
beta.natural.to.params = function(x) {
  p = x + 1
  dimnames(p)[[1]] = c("a", "b")
  p
}

# Converts from beta parameters to moments.
beta.params.to.moments = function(p) {
  p1 = array(p, dim=c(2, length(a)/2))
  a = p1[1,]
  b = p1[2,]

  m = array(cbind(   ,    ), dim=dim(a))
  dimnames(m) = dimnames(a)
  dimnames(m)[[1]] = c("m", "v")
  m
}

# Converts from beta moments to parameters.
beta.moments.to.params = function(m) {




}

# Distribution of x ~ N(m,v) | Ax ~ N(b,b.var).
lin.constraint.1 = function(m, v, A, b, b.var) {

  M = A %*% (v * t(A))
  diag(M) = diag(M) + b.var
# cat("M =", M, "\n")

  # using Cholesky would be slightly faster, but maybe less stable
  M.qr = qr(M)

  r = cbind(m = m - v * as.vector(t(A) %*% solve(M.qr, A %*% m - b)),
    v = v - v * apply(A * solve(M.qr, A), 2, sum) * v)
# print(colnames(r))
  r
}

# Posterior for several independent linear constraints.
# Args:
#   m, v - mean and variance
#   A, b - these give individual constraints
# Returns: list with elements "m" and "v" (the posterior)
lin.constraint.1eq = function(m, v, A, b) {

  # hopefully none of these end up <= 0...
  M = apply(A * v * A, 1, sum)
# cat("M =", M, "\n")

  z = (apply(A * m, 1, sum) - b) / M
# print(z)
  y = A * (A / M)
# print(y)

  list(m = t( t(m) - t(v * A * z) ),
    v = v - v * y * v )
}

# Estimates messages from each factor.
# Args:
#   m - the incoming messages
#   a, b - these give the constraint that ax = b.
# Returns: the outgoing messages
approx.region.beta.message = function(m, a, b) {



}

# Approximates a region using beta distributions.
# Args:
#   a, b - these give the constraint that ax = b.
#     All entries of each should be in [0,1] .
# Returns: list with elements
#   m, v - mean and variance of the posterior estimate
#   x - parameters of the posterior ("a" and "b")
approx.region.beta = function(a, b, max.iters = 20) {


  # messages from each factor
  q1 = array(0, dim=c(2, nrow(a), ncol(a)))
  dimnames(q1)[[1]] = c("e1", "e2")

  # posterior
  q = apply(q1, 2, sum)

  for(iter in 1:max.iters) {



    # update messages from each factor
#    q1 = q1 + 


    # update posterior
    q = apply(q1, 2, sum)
  }


}








