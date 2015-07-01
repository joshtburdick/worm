# Utilities for the gamma distribution.

source("git/unmix/ept/beta.r")

# Utility to allow writing conversion functions,
# with the parameters as the first dimension.
# ??? should the "f" functions be rewritten to consider
# the rownames to be the parameter names?
array.apply = function(f) function(a) {
  x1 = matrix(a, nrow=dim(a)[[1]])
  rownames(x1) = dimnames(a)[[1]]
  r = t(f(t(x1)))
  a = array(r, dim=dim(a), dimnames=dimnames(a))
  dimnames(a)[[1]] = rownames(r)
  a
}

# Converts from standard to natural parameters.
gamma.s2n = array.apply(function(a) {
  cbind(e1 = a[,"a"] - 1, e2 = - a[,"b"])
})

# Converts from natural to standard parameters.
gamma.n2s = array.apply(function(a) {
  cbind(a = a[,"e1"] + 1, b = - a[,"e2"])
})

# Converts from mean and variance to standard parameters (shape and rate.)
gamma.mv2s = array.apply(function(a) {
  b = a[,"m"] / a[,"v"]
  cbind(a = a[,"m"] * b, b = b)
})

# Converts from standard parameters (shape and rate) to mean and variance.
gamma.s2mv = array.apply(function(a) {
  r = cbind(m = a[,"a"] / a[,"b"],
    v = a[,"a"] / (a[,"b"]^2))
  r
})

# More utilities.
gamma.mv2n = function(a) gamma.s2n(gamma.mv2s(a))
gamma.n2mv = function(a) gamma.s2mv(gamma.n2s(a))

# Marginal mean and variance of a Dirichlet.
dirichlet.s2mv = function(a) {
  a0 = sum(a)
  beta.s2mv(rbind(a = a, b = a0 - a))
}

# Marginal mean and variance of independent gamma variates,
# conditional on them summing to 1.
# (Not sure this is right.)
# Args:
#   x - gamma distribution (natural parameters)
# Returns: natural parameters of gamma, conditional on that.
gamma.cond.1 = function(x) {

  # convert to standard parameters
  xs = gamma.n2s(x)

  # compute Dirichlet moments
  a = dirichlet.s2mv(xs["a",])

  # scale by the rate parameter of the original gamma variables
  r = 1 / xs["b",]
  a1 = rbind(m = r * a["m",], v = r * r * a["v",])

  # normalize
  s = 1 / sum(a1["m",])
  gamma.mv2n(rbind(m = s * a1["m",], v = s * s * a1["v",]))
}

# Marginal mean and variance of independent gamma variates,
# conditional on a linear constraint.
# Args:
#   x - gamma distribution (natural parameters)
#   a, b - these give the constraint that ax = b
#     (we assume a, b > 0)
# Returns: natural parameters of gamma, conditional on that.
gamma.cond = function(x, a, b) {

  xs = gamma.n2s(x)
  xs["b",] = xs["b",] * a

  y = gamma.n2s(gamma.cond.1(gamma.s2n(xs)))

  y["b",] = y["b",] / b
  gamma.s2n(y)
}

#   Estimates marginals of independent gamma variates,
# conditional on them summing to 1.
#   This is an approximation; but if the gammas all have
# the same scale, then it amounts to moment-matching
# a beta with a gamma (which should often be good.)
# Args:
#   x - gamma distributions (as natural parameters)
# Returns: same, conditional on them summing to 1
# XXX if the scales are different, then the result of
# this may not sum to 1, which seems incorrect
gamma.conditional.approx.1.beta = function(x) {

  xa = gamma.n2s(x)

  # compute sum of all but one of these
  x.mv = gamma.n2mv(x)
#  x.mv[ is.nan(x.mv) ] = 0    # ignore NaN and NA
  xb = gamma.mv2s( apply(x.mv, 1, sum) - x.mv )
#  xb[ is.nan(xb) ] = 0    # ignore NaN and NA

  # if we pretend that xa and xb have the same rate,
  # then xa and xb now define a beta distribution (with
  # parameters given by their shapes); moment-match
  # this as a gamma
  b = rbind(a=xa["a",], b=xb["a",])
  mm = gamma.mv2s(beta.s2mv(b))

  # lastly, xa and xb, in general, will have different
  # rates; correct for this difference (essentially
  # multiplying by an exponential distribution)
  mm["b",] = mm["b",] + xa["b",] - xb["b",]

  gamma.s2n(mm)
}

#   Estimates marginals of independent gamma variables,
# conditional on them summing to 1. This is a simpler
# method than gamma.conditional.approx.1.beta .
# XXX and incorrect, I think.
# Args:
#   x - gamma distributions (as natural parameters)
# Returns: same, conditional on them summing to 1.
gamma.conditional.approx.1.rescale = function(x) {
  # XXX this only works on one at a time

  # compute sum of these
  x.mv = gamma.n2mv(x)
  s = sum(x.mv["m",])




  list(x.mv = x.mv)
}

# Approximates marginals of gamma-distributed variables,
# given constraints on their expected sums.
# Args:
#   x - natural parameters of gamma distributions
#   A, b - these give the constraints that A * x = b
# Returns: x, but conditional on the constraint
# Deprecated; I don't think this is correct.
gamma.conditional.approx.array = function(x, A, b) {

  # determine how much to scale each
  s1 = b / A
  s1[A==0] = 0
# print(s1)

  # scale these
  xs = x
  xs["e2",,] = xs["e2",,] * t(s1)

  # zero out cases in which A = 0
  xs["e1",,][ t(A==0) ] = 0
  xs["e2",,][ t(A==0) ] = -Inf

# print(xs["e2",,])
# print(gamma.n2mv(xs))

  # condition on the scaled distribution summing to 1
  xsc = array(apply(xs, c(3), gamma.conditional.approx.1.beta),
    dim=dim(xs), dimnames=dimnames(x))
# print(xsc)
  dimnames(xsc) = dimnames(x)

  # undo the scaling
  xc = xsc
  xc["e2",,] = xsc["e2",,] / t(s1)

  # cases in which A = 0 aren't affected by the constraint
  xc["e1",,][ t(A==0) ] = x["e1",,][ t(A==0) ]
  xc["e2",,][ t(A==0) ] = x["e2",,][ t(A==0) ]

  xc
}


# The density function for the gamma, conditional on vars summing to 1.
# XXX not sure this is right.
# ETA: I think it's correct for two variables, but not otherwise.
gamma.conditional.density = function(p) function(x) {
  r = dgamma(x, shape=p["a",1], rate=p["b",1])
  for(i in 2:ncol(p)) {
    r = r * dgamma((1-x), shape=p["a",i], rate=p["b",i])
  }
  r
}

# Another attempt at the gamma conditional thing.
# Returns: mean and variance of first variable, conditional
# on all of them summing to 1.
gamma.conditional.numerical.1 = function(p) {
  f = gamma.conditional.density(p)
rull
  # compute the moments
  # XXX assuming successful integration
  x0 = integrate(f, 0, 1)$value
  x1 = integrate(function(x) x * f(x), 0, 1)$value
  x2 = integrate(function(x) x * x * f(x), 0, 1)$value

  rbind(m = x1 / x0, v = x2/x0 - (x1/x0)^2)
}

# (Hopefully) convenient interface for the above.
gamma.conditional.numerical = function(a)
  gamma.mv2n(gamma.conditional.numerical.1(gamma.n2s(a)))




