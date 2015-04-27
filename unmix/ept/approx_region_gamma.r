# Approximates a region using gamma distributions.

library("limSolve")

source("git/unmix/ept/beta.r")
source("git/unmix/ept/gamma.r")

# Finds the m'th moment of a gamma truncated to
# some interval. Possibly not used.
# Args:
#   m - the moment to find
#   x - where to truncate the distribution (should be >= 0)
#   a, b - parameters of gammas
# Returns: the m'th moment
trunc.gamma.moment = function(m, x, a, b)
  exp( lgamma(a+m) - lgamma(a) - m*log(b) ) * pgamma(x, a+m, b)

# Converts from mean and variance to moments, and back.
mv2m = function(x) {
  r = x
  dimnames(r)[[1]] = c("x1", "x2")
  r["x2",,] = x["v",,] + x["m",,]^2
  r
}
m2mv = function(x) {
  r = x
  dimnames(r)[[1]] = c("m", "v")
  r["v",,] = x["x2",,] - x["x1",,]^2
  r
}

# Incorporates some linear constraints, using a somewhat
# vague and experimental moment-matching method, which may
# well be completely incorrect.
# XXX ... and indeed, it doesn't seem to work.
# Args:
#   x - gamma distributions (as natural parameters)
#   A, b - these give constraints on the sums of x
# Returns: moment-matched estimate
gamma.conditional.approx.mm = function(x, A, b) {
  # number of constraints
  n = length(b)

  # compute the posteriors for each term
  x1 = gamma.conditional.approx(x, A, b)

  # incorporate constraint, using a weighted combination of moments
#  x2.moments = (mv2m(gamma.n2mv(x)) + mv2m(gamma.n2mv(x1))) / 2
  x2.moments = ((n-1) * mv2m(gamma.n2mv(x)) + mv2m(gamma.n2mv(x1))) / 2

  # return that, converted to natural parameters
  gamma.mv2n(m2mv(x2.moments))
}

# Approximates a region subject to Ax = b, x >= 0.
#   A, b - these define exact constraints.
# Returns: list with components
#   x - prediction mean and variance
approx.region.gamma = function(A, b, max.iters=100) {

  # the prior (initially something non-flat.)
  x1 = rbind(m = rep(1, ncol(A)), v = rep(1^2, ncol(A)))
  prior = gamma.mv2n(x1)
  colnames(prior) = colnames(A)
prior = prior / 2

  # the approximating terms
  term = array(0, dim=c(2, ncol(A), nrow(A)),
    dimnames = list(c("e1", "e2"), colnames(A), rownames(A)))
  # initialize with all terms having mean = 1, variance = 1
  # ??? or flat?
  term["e1",,] = 0
  term["e2",,] = -1

  # posterior is the sum of terms
  # FIXME start with truncated pseudoinverse?
  q = apply(term, c(1,2), sum) + prior

  for (iter in 1:max.iters) {

    # sum of messages from all other terms
    q1 = - sweep(term, c(1,2), q, "-")
# print(gamma.n2mv(q1))

    # messages from each factor, conditional on constraint
    q.new = gamma.conditional.approx(q1, A, b)

    # difference between those messages, and current posterior
    change = sweep(q.new, c(1,2), q, "-")
#    change = q1 - q.new

    # update messages from each factor (possibly damped)
    term = term + 1 * change

    # update posterior, as sum of messages
    q = apply(term, c(1,2), sum) + prior

    # FIXME check for convergence

  }

  list(x = q, term = term)
}

# Another attempt at the same thing.
approx.region.gamma.2 = function(A, b, max.iters=100) {

  # the prior (initially something non-flat.)
  x1 = rbind(m = rep(1, ncol(A)), v = rep(1, ncol(A)))
  prior = gamma.mv2n(x1)
  colnames(prior) = colnames(A)

  # the approximating terms
  term = array(0, dim=c(2, ncol(A), nrow(A)),
    dimnames = list(c("e1", "e2"), colnames(A), rownames(A)))
  # initialize with mean = 1, variance = 1
  term["e1",,] = 0
  term["e2",,] = -1

  # posterior is the sum of terms
  # FIXME start with truncated pseudoinverse?
  q = prior + apply(term, c(1,2), sum)

  for (iter in 1:max.iters) {

    # sum of messages from all other terms
    q1 = - sweep(term, c(1,2), q, "-")
# print(gamma.n2mv(q1))

    # messages from each factor, conditional on constraint
    q.new = gamma.conditional.approx(q1, A, b)

    # difference between those messages, and current posterior
    change = - sweep(q.new, c(1,2), q, "-")
#    change = q1 - q.new

    # update messages from each factor (possibly damped)
    term = term + 1 * change

    # update posterior, as sum of messages
    q = prior + apply(term, c(1,2), sum)

    # FIXME check for convergence
print(q)
  }

  list(x = q, term = term)
}

# Wrapper for gamma-conditional-approx which avoids the
# whole multidimensional array thing.
gamma.conditional.approx.vec = function(x, A, b) {

}

if (TRUE) {

A0 = t(rep(1,4))
r = approx.region.gamma(A0, 1, max.iters=100)
# this exactly matches the moments of the actual marginals,
# which are ~ Beta(1,4) (for what that's worth)
# to see this, eval
#     gamma.n2mv(r$x)

x0 = array(gamma.s2n(rbind(a=c(2,3,4,5), b=rep(1,4))), dim=c(2,4,1))
dimnames(x0)[[1]] = c("e1", "e2")
r.gca = gamma.conditional.approx(x0, A0, 1)
# this seems fairly plausible

A1 = rbind(c(1,1,0), c(0,1,1))
r1 = approx.region.gamma(A1, c(1,1), max.iters = 200)

A2 = rbind(
  c(1,1,1,1,1,0,0),
  c(0,0,0,1,1,1,1))
rownames(A2) = c("ceh-6", "hlh-1")
colnames(A2) = c("ABal", "ABar", "ABpl", "ABpr", "MS", "E", "P2")   # XXX bogus
r2 = approx.region.gamma(A2, c(0.5,0.5), max.iters=100)

}

