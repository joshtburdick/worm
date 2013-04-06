# Approximates a region using gamma distributions.

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


# Approximates a region subject to Ax = b, x >= 0.
#   A, b - these define an exact constraint.
# Returns: list with components
#   x - prediction mean and variance
approx.region.gamma = function(A, b, max.iters=100) {

  # the prior (initially something non-flat.)
  x1 = rbind(m = rep(1, ncol(A)), v = rep(1, ncol(A)))
  prior = gamma.mv2n(x1)
  colnames(prior) = colnames(A)
prior = prior / 1

  # the approximating terms
  term = array(0, dim=c(2, ncol(A), nrow(A)),
    dimnames = list(c("e1", "e2"), colnames(A), rownames(A)))
  term["e1",,] = 0
  term["e2",,] = -1

  # posterior is the sum of terms
  # FIXME start with truncated pseudoinverse?
  q = prior + apply(term, c(1,2), sum)

  for (iter in 1:max.iters) {

    # sum of messages from all other terms
    q1 = - sweep(term, c(1,2), q, "-")
# print(gamma.n2mv(m.to))

    # messages from each factor, conditional on constraint
    q.new = gamma.conditional.approx(q1, A, b)

    # difference between those messages, and current posterior
#    change = sweep(m, c(1,2), q, "-")

    # update messages from each factor (possibly damped)
#    term = term + 1 * change
    term = q.new - q1;

    # update posterior, as sum of messages
    q = prior + apply(term, c(1,2), sum)

    # FIXME check for convergence

  }

  list(x = q, term = term)
}

A0 = t(rep(1,4))
x0 = array(gamma.s2n(rbind(a=c(2,3,4,5), b=rep(1,4))), dim=c(2,4,1))
dimnames(x0)[[1]] = c("e1", "e2")
# r = gamma.conditional.approx(x0, A0, 1)

r = approx.region.gamma(A0, 1, max.iters=100)

A1 = rbind(c(1,1,0), c(0,1,1))
r1 = approx.region.gamma(A1, c(1,1), max.iters = 6)

if (FALSE) {
A2 = rbind(
  c(1,1,1,1,1,0),
  c(0,0,0,1,1,1))
rownames(A) = c("ceh-6", "hlh-1")
colnames(A) = c("ABal", "ABar", "ABpl", "ABpr", "EMS", "P2")   # XXX bogus
r2 = approx.region.gamma(A2, c(1,1), max.iters=3)
}

