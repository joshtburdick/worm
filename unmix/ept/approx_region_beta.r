# Attempts to approximate a region using beta distributions
# (instead of Gaussians.)



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

# Approximates a region using beta distributions.
# Args:
#   a, b - these give the constraint that a x = b.
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








  }


}








