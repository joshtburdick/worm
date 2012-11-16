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
  stopifnot(all(dimnames(p)[[1]] == c("a", "b")))
  f = function(x) {
    a = x[1]
    b = x[2]
    c(mu = a / (a+b), si2 = a*b / ((a+b)^2 * (a+b+1)))
  }
  m = apply(p, c(-1), f)
  dimnames(m)[[1]] = c("mu", "si2")
  m
}

# Converts beta moments to parameters.
# from http://en.wikipedia.org/wiki/Beta_distribution
beta.moments.to.params = function(m) {
  stopifnot(all(dimnames(m)[[1]] == c("mu", "si2")))
  f = function(x) {
    mu = x[1]
    si2 = x[2]
    s = (mu * (1-mu) / si2) - 1
    c(a = mu * s, b = (1-mu) * s)
  }
  p = apply(m, c(-1), f)
  dimnames(p)[[1]] = c("a", "b")
  p
}

# Normal moment matching (array version.)
normal.moment.match.a = function(p) {
  stopifnot( length(dim(p)) == 3 )
  stopifnot(all(dimnames(p)[[1]] == c("mu", "si2")))
  m = p["mu",,]
  v = p["si2",,]

  m1 = -as.vector(m)
# v[v<0] = 1e10   # XXX hack
  s = sqrt(v)
  z = -m1 / s
  a = dnorm(z) / pnorm(z)

  # hack to deal with when z is very negative
  r = c(m = ifelse(z < -30, 0, - (m1 - s * a)),
    v = ifelse(z < -30, 0, v * (1 - z*a - a^2)))
#  r = c(m = - (m1 - s * a), v = v * (1 - z*a - a^2))

  p["mu",,] = r$m
  p["si2",,] = r$v
  p
}

# Normal moment matching (slow, flexible version.)
normal.moment.match.slow = function(p) {
  stopifnot(all(dimnames(p)[[1]] == c("mu", "si2")))

  f = function(x) {
    m = x["mu"]
    v = x["si2"]
    m1 = -as.vector(m)
    # v[v<0] = 1e10   # XXX hack
    s = sqrt(v)
    z = -m1 / s
    a = dnorm(z) / pnorm(z)

    # hack to deal with when z is very negative
    c(m = ifelse(z < -30, 0, - (m1 - s * a)),
      v = ifelse(z < -30, 0, v * (1 - z*a - a^2)))
    #  c(m = - (m1 - s * a), v = v * (1 - z*a - a^2))
  }

  p = apply(p, -1, f)
  dimnames(p)[[1]] = c("mu", "si2")
  p
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
# ("what each factor would like the answer to be")
approx.region.beta.message = function(m, a, b) {



}

# Approximates a region using beta distributions.
# Args:
#   a, b - these give the constraint that ax = b.
#     All entries of each should be in [0,1] .
# Returns: list with elements
#   m, v - mean and variance of the posterior estimate
#   x - parameters of the posterior ("a" and "b")
approx.region.beta.1 = function(a, b, max.iters = 5) {

  # messages from each factor
  q1 = array(0, dim=c(2, nrow(a), ncol(a)))
  dimnames(q1)[[1]] = c("e1", "e2")
  # ??? is this a good initial setting?

  # posterior
print(list(q1=q1))
  q = apply(q1, c(1,3), sum)
print(list(q=q))
  for(iter in 1:max.iters) {

    # compute incoming messages to each factor
    q1.in = q1
    for(j in 1:(dim(q1.in)[2]))
      q1.in[,j,] = q - q1[,j,]

    # compute outgoing messages from each factor
    q1.mv = beta.params.to.moments(beta.natural.to.params(q1.in))
print(list(q1.mv=q1.mv))
    r = lin.constraint.1eq(t(q1.mv["mu",,]), t(q1.mv["si2",,]), a, b)
print(list(r=r))
    q2 = q1.mv
    q2["mu",,] = t(r$m)
    q2["si2",,] = t(r$v)

    # moment-match each of these
    q.mm = beta.params.to.natural(beta.moments.to.params(normal.moment.match.slow(q2)))
print(list(q.mm=q.mm))
    # update messages from each factor
#    q1 = q1 + (q.mm - q)
    # update posterior
    for(j in 1:(dim(q1)[2]))
      q1[,j,] = q1[,j,] + (q.mm[,j,] - q)
#      q1[,j,] = q1[,j,] + (q.mm[,j,] - q)

    q = apply(q1, 2, sum)
  }

  list(q1 = q1, q = q)
}

# toy test problems
r1 = approx.region.beta.1(t(c(1,1,1)), 1)

