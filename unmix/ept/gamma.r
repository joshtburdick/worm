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

# Finds the marginal distributions of independant gamma variables,
# conditional on their sum adding up to exactly one.
# (Analogous to lin.constraint() in the EP code.)
# ??? is this correct? Presumably their sum should be
# just the sum of the expected values?
# Currently only handles a single constraint.
# Args:
#   x - natural parameters of gamma distributions
#   (should be an array with three dimensions,
#   first being "e1" and "e2")
# Returns: x, conditional on the x's summing to 1
gamma.conditional.1 = function(x) {

  # first, compute sum of all x's (by scaling them
  # all "blurred" together) ??? is this correct?
#  s = gamma.n2mv(as.matrix(apply(x, 1, sum)))
# * c(ncol(x), ncol(x)^2)
# print(s)
  # compute mean
#  m = s["m",1] * ncol(x)

  x.mv = gamma.n2mv(x)

  # compute mean of each
  m = apply(x.mv["m",,], 1, sum)
print(m)

  # scale by that mean
  m1 = rbind(m, m*m)
  r = sweep(x.mv, c(1,3), as.vector(m1), "/")

  r = gamma.n2mv(x) / c(m, m*m)
print(r)
  # return that, moment-matched as gamma
  gamma.mv2n(r)
}

# Finds the marginal distributions of independant gamma variables,
# conditional on their sum adding up to exactly one.
# (Analogous to lin.constraint() in the EP code.)
# ??? is this correct?
# Args:
#   x - natural parameters of gamma distributions
# Returns: x, conditional on the x's summing to 1
gamma.conditional.beta.1 = function(x) {

  # first, moment-match these as beta distributions
  xb = beta.mv2n(gamma.n2mv(x))
print(round(xb,3))

  # then, compute one minus each of these
  # (this is the same as switching the a and b parameters)
  xb1 = rbind(e1 = xb["e2",], e2 = xb["e1",])
print(round(xb1,3))

  # compute average of all but one of these
print(ncol(x))
  r = (apply(xb1, 1, sum) - xb1) # / (ncol(x)-1)
print(round(r, 3))

  # re-swap
#  r1 = rbind(e1 = r["e2",], e2 = r["e1",])
#print(round(r1, 3))

  # average this with original marginal
  p = (r + xb)   # / 2
print(round(p, 3))

  # return that, moment-matched as gamma, and rescaled
  gamma.mv2n(beta.n2mv(p))
}

# Variation of this which breaks the product into an
# exponential distribution times a beta, moment-matches
# the beta as a gamma, then multiplies by the
# exponential distribution.
# Args:
#   x - gamma distributions (as natural parameters)
# Returns: same, conditional on them summing to 1
gamma.conditional.beta.2 = function(x) {

  xa = gamma.n2s(x)

  # compute product of all but one of these
  xb = gamma.n2s( apply(x, 1, sum) - x )

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

# As above, but allows arbitrary linear constraints.
# Args:
#   x - natural parameters of gamma distributions
#   A, b - the constraint that Ax = b.
#     A must be non-negative, and b must be positive.
# Returns: parameters x, conditional on constraint.
gamma.conditional = function(x, A, b) {
  # ignore cases in which b = 0
  x1 = x[ , b > 0 ]
  A1 = A[ b > 0 , ]
  b1 = b[ b > 0 ]

  # determine how much to scale each
  s1 = A1 / b1

  # scale (in mean-and-variance space)
  x1.scaled = gamma.mv2n( t( t(gamma.n2mv(x1)) * c(s1, s1*s1) ) )  

  # condition on x1 summing to 1
  x1.c = gamma.conditional.1(x1.scaled)

  # undo the scaling
  x1.c.unscaled = gamma.mv2n( t( t(gamma.n2mv(x1.c)) / c(s1, s1*s1) ) )

  # store cases where b > 0
  x.c = x
  x.c[ , b > 0 ] = x.c1
  x.c
}

# The density function for the gamma, conditional on vars summing to 1.
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
gamma.conditional.numerical = function(p) {
  f = gamma.conditional.density(p)

  # compute the moments
  # XXX assuming successful integration
  x0 = integrate(f, 0, 1)$value
  x1 = integrate(function(x) x * f(x), 0, 1)$value
  x2 = integrate(function(x) x * x * f(x), 0, 1)$value

  rbind(m = x1 / x0, v = x2/x0 - (x1/x0)^2)
}

# various tests for gamma.conditional
x1.mv = beta.s2mv(rbind(a=c(1,1,1), b=c(1,2,3)))
x1 = gamma.mv2n(x1.mv)

x2 = gamma.mv2n(rbind(m=c(0.4,0.2,0.3), v=c(0.08,0.02,0.08)))

x3.mv = array(c(2,1, 2,1, 2,2,
    3,3, 3,1, 1,1,
    1,1, 1,1, 1,1),
  dim=c(2,3,3),
  dimnames=list(c("m","v"), c("x1","x2","x3"), c("a","b","c")))
x3 = gamma.mv2n(x3.mv)


