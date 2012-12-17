# Attempts to approximate a region using beta distributions
# (instead of Gaussians.) It's a fairly conservative design,
# in that other than that, it's basically the same.

# Possibly deprecated (Dirichlet seems like the thing to
# use eventually, anyway.)

source("git/unmix/ept/beta.r")



# Computes distribution of a beta.
# Args:
#   x - the distribution (as a set of beta distributions.)
#   a, b - these give the constraint that ax = b.
# Returns: the distribution of x, conditional on ax = b.
beta.conditional.on.sum = function(x, a, b) {





}

# Approximates a region using a beta distribution
# for each cell.
#   A, b - these give the constraint that Ax = b.
# Returns: list with components
#   x - prediction mean and variance
approx.region.beta = function(A, b, max.iters=100) {

  # the posterior (initially something non-flat.)
  # FIXME start with truncated pseudoinverse?
  x1 = rbind(m = rep(1, ncol(A)), v = rep(1, ncol(A)))
  x = gamma.mv2n(x1)
  colnames(x) = colnames(A)

  # the messages from each factor (initially flat)
  m = array(0, dim=c(2, ncol(A), nrow(A)),
    dimnames = list(c("e1", "e2"), colnames(A), rownames(A)))

  for (iter in 1:max.iters) {

    # sum of messages from all other factors
    m.to = - sweep(m, c(1,2), x, "-")
print(gamma.n2mv(m.to))

    # messages from each factor, conditional on constraint
    m.from = gamma.conditional(m.to, A, b)

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





# toy test problems
r1 = approx.region.beta.1(t(c(1,1,1)), 1)

