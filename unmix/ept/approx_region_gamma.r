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

  # the posterior (initially something non-flat.)
  # FIXME start with truncated pseudoinverse?
  x1 = rbind(m = rep(1, ncol(A)), v = rep(1, ncol(A)))
  x = gamma.mv2n(x1)
  colnames(x) = colnames(A)

  # the approximating terms
  m = array(0, dim=c(2, ncol(A), nrow(A)),
    dimnames = list(c("e1", "e2"), colnames(A), rownames(A)))
  m["e1",,] = 0
  m["e2",,] = -1

  # posterior is the sum of terms
  x = apply(m, c(1,2), sum)

  for (iter in 1:max.iters) {

    # sum of messages from all other factors
    m.to = - sweep(m, c(1,2), x, "-")
print(gamma.n2mv(m.to))

    # messages from each factor, conditional on constraint
    m.from = gamma.conditional.approx(m.to, A, b)

    # difference between those messages, and current posterior
    m.change = sweep(m.from, c(1,2), x, "-")

    # update messages from each factor (possibly damped)
    m = m + 1 * m.change

    # update posterior, as sum of messages
    x = apply(m, c(1,2), sum)

    # FIXME check for convergence

  }

  list(x = x, m = m)
}


A = rbind(
  c(1,1,1,1,1,0),
  c(0,0,0,1,1,1))
rownames(A) = c("ceh-6", "hlh-1")
colnames(A) = c("ABal", "ABar", "ABpl", "ABpr", "EMS", "P2")   # XXX bogus

A0 = t(rep(1,5))
x0 = array(gamma.s2n(rbind(a=rep(1,5), b=rep(1,5))), dim=c(2,5,1))
dimnames(x0)[[1]] = c("e1", "e2")
# r = gamma.conditional(x0, A0, 1)



# r = approx.region.gamma(A0, c(1), max.iters=10)


