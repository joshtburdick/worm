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
# ??? is this correct?
# Args:
#   x - natural parameters of gamma distributions
# Returns: x, conditional on the x's summing to 1
gamma.conditional.1 = function(x) {

  # first, moment-match these as beta distributions
  xb = beta.mv2n(gamma.n2mv(x))
print(round(xb,3))

  # then, compute one minus each of these
  # (this is the same as switching the a and b parameters)
  xb1 = rbind(e1 = xb["e2",], e2 = xb["e1",])
print(round(xb1,3))

  # compute average of all but one of these
print(ncol(x))
  r = (apply(xb1, 1, sum) - xb1) / (ncol(x)-1)
print(round(r, 3))

  # re-swap
  r1 = rbind(e1 = r["e2",], e2 = r["e1",])
print(round(r1, 3))

  # average this with original marginal
  p = (r1 + xb) / 2
print(round(p, 3))

  # return that, moment-matched as gamma, and rescaled
  gamma.mv2n(beta.n2mv(p))
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

x1.mv = beta.s2mv(rbind(a=c(1,1,1), b=c(1,1,1)))
x1 = gamma.mv2n(x1.mv)

x2 = gamma.mv2n(rbind(m=c(0.5,0.5), v=c(0.1,0.1)))

