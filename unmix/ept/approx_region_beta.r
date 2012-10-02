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

# Converts beta parameters to moments.
# from http://en.wikipedia.org/wiki/Beta_distribution
beta.params.to.moments = function(p) {
  a = p["a",]
  b = p["b",]
  m = rbind(mu = a / (a+b), si2 = a*b / ((a+b)^2 * (a+b+1)) )
#  rownames(m) = c("mu", "si2")
  m
}

# Converts beta moments to parameters.
# from http://en.wikipedia.org/wiki/Beta_distribution
beta.moments.to.params = function(m) {
  mu = m["mu",]
  si2 = m["si2",]
  s = (mu * (1-mu) / si2) - 1     # ??? not sure what this is
  m = rbind(a = mu * s, b = (1-mu) * s)
#  rownames(m) = c("a", "b")
  m
}


# Converts from beta parameters to moments.
# FIXME this would preserve array dimensions, which
# might be nice. But may be irrelevant.
#beta.params.to.moments = function(p) {
#  p1 = array(p, dim=c(2, length(a)/2))
#  a = p1[1,]
#  b = p1[2,]

#  m = array(cbind(   ,    ), dim=dim(a))
#  dimnames(m) = dimnames(a)
#  dimnames(m)[[1]] = c("m", "v")
#  m
#}

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








