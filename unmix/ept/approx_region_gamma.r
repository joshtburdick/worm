# Approximates a region using gamma distributions.

source("git/unmix/ept/gamma.r")

# Posterior for several independent linear constraints.
# Args:
#   m, v - mean and variance
#   A, b - these give individual constraints
# Returns: list with elements "m" and "v" (the posterior)
lin.constraint.1eq = function(m, v, A, b) {

  # hopefully none of these end up <= 0 (this should
  # be the case if at least one entry in each group
  # of cells has nonzero variance)
  M = apply(A * v * A, 1, sum)

# cat("M =", M, "\n")

  z = (apply(A * m, 1, sum) - b) / M
# print(z)
  y = A * (A / M)
# print(y)

  list(m = t( t(m) - t(v * A * z) ),
    v = v - v * y * v )
}

# Does unmixing.
#   A, b - these define an exact constraint.
# Returns: list with components
#   x - prediction mean and variance
approx.region.gamma = function(A, b, max.iters=100) {
  
  # the messages from each factor
  # (initially epsilon. FIXME start with truncated pseudoinverse?)
  m = array(0, dim=c(2, ncol(A), nrow(A)),
    dimnames = list(c("e1", "e2"), colnames(A), rownames(A)))

  # initialize moment-matched estimate (in terms of mean and variance)
  mm = 0 * m
  dimnames(mm)[[1]] = c("m", "v")

  # the posterior (initially something non-flat.)
  # FIXME start with truncated pseudoinverse?
  x1 = rbind(m = rep(1, ncol(A)), v = rep(1, ncol(A)))
  x = gamma.mv2n(x1)
  colnames(x) = colnames(A)

  for (iter in 1:max.iters) {

    # compute messages to each factor (in terms of "mean" and "variance")
    m.to = gamma.n2mv(as.vector(x) - m)
print(m.to)

    # compute posterior "from" each factor
    mv = lin.constraint.1eq(t(m.to["m",,]), t(m.to["v",,]), A, b)
print(mv)

    # normal moment-match

    mm["m",,] = t(mv$m)
    mm["v",,] = t(mv$v)
print(mm)

    # update messages from each factor (damped)
    m.change = 0.5 * ( gamma.mv2n(mm) - as.vector(x) )
    m = m + m.change

    # update posterior, as sum of messages
    x = apply(m, c(1,2), sum)
  }

  list(x = x, m = m)
}


A = rbind(
  c(1,1,1,1,1,0),
  c(0,0,0,1,1,1))
rownames(A) = c("ceh-6", "hlh-1")
colnames(A) = c("ABal", "ABar", "ABpl", "ABpr", "EMS", "P2")   # XXX bogus




