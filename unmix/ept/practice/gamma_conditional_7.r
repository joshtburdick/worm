# Attempt at finding marginals of
#   x ~ Gamma(a_i, b_i) | Ax = b
# this time by rescaling the (presumed) corners
# of Ax=b to a unit simplex.

library(corpcor)
library(MASS)

source("git/unmix/ept/gamma.r")
source("git/unmix/ept/practice/gamma_conditional_6.r")
source("git/unmix/ept/practice/simplex_corner.r")

# One approximation of marginals of
#   x ~ Gamma(s_i, r_i) | sum(x) = 1
# by basically rescaling all the parameters.
# XXX trying gamma.conditional.approx.1.beta instead.
# Args:
#   x - gamma parameters of x (as natural parameters)
# Returns: natural parameters of x | sum(x) = 1
gamma.cond.sum1 = function(x) {

  # convert to mean and variance
  x1 = gamma.n2mv(x)
  s = sum(x1["m",])

  # XXX there's probably a better way to do this
  x1s = x1
  x1s["m",] = x1s["m",] / s
  x1s["v",] = x1s["v",] / (s*s)

  # convert back to natural parameters
  gamma.mv2n(x1s)
}

# Approximates marginals of
#   x ~ Gamma(s_i, r_i) | Ax = b
# by rescaling to the above form, solving, then
# scaling back.
# Args:
#   A, b - these give the linear system
#   x - gamma parameters of x (as natural parameters)
# Returns: natural parameters of x | sum(x) = 1
gamma.cond.1 = function(A, b, x) {

  x.mv = gamma.n2mv(x)

  S = simplex.corner(A, b)
  S.inv = pseudoinverse(S)

  m.s = as.vector(S.inv %*% x.mv["m",])
  # XXX trying just using the diagonal of the covariance,
  # which may be an oversimplification
  # XXX ignoring memory usage
  v.s = diag( S.inv %*% diag(x.mv["v",]) %*% t(S.inv) )

  # marginals of the transformed problem
  # (was using gamma.cond.sum1 earlier)
  c1 = gamma.n2mv(gamma.conditional.approx.1.beta(gamma.mv2n(
    rbind(m = m.s, v = v.s))))

  m.p = as.vector(S %*% (c1["m",]))
  v.p = diag( S %*% diag(c1["v",]) %*% t(S) )

  gamma.mv2n(rbind(m = m.p, v = v.p))
}

A = matrix(rgamma(10, shape=1, rate=1), nrow=2)

# Test of the unit simplex thing.
# Args:
#   x - gamma variables (natural parameters)
# Returns: x | sum(x) = 1, as a list with two elements
#   as mean and variance, x.predicted and x.sampled
test.unit.simplex = function(x) {
  list(x.predicted = gamma.n2mv(gamma.cond.sum1(x)),
    x.sampled = gamma.n2mv(gamma.cond.sampling(
      x, t(rep(1, ncol(x))), 1)))
}

x1 = gamma.s2n(rbind(a=c(1,1,1), b=c(1,1,1)))
r1 = gamma.cond.1(t(c(1,2,3)), 1, x1)

A3 = rbind(c(1,1,0), c(0,1,1))
r3 = gamma.cond.1(A3, c(1,2), x1)

