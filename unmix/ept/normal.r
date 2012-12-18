# Normally-distributed variables, and related utilities.

# Creates a multivariate normal vector variable.
norm.var = function(n) {
  b = matrix(0, nrow=n, ncol=2)
  # XXX I didn't know what to call these, but named them this, since
  # Wikipedia calls the canonical parameters "eta"
  # ??? is "n" needed?
  colnames(b) = c("e1", "e2")
  list(b = b, observed=FALSE)
}

# Converts from mean and variance to canonical parameters.
mean.and.variance.to.canonical = function(a) {
  r = cbind(e1 = a[,"m"] / a[,"v"], e2 = -0.5 / a[,"v"])
  # deal with special case when m = v = 0.
  r[ (a[,"m"]==0 & a[,"v"]==0), "e1" ] = 0
  r[ (a[,"m"]==0 & a[,"v"]==0), "e2" ] = -Inf
  r
}

# Converts from canonical parameterization to mean and variance.
canonical.to.mean.and.variance = function(a) {
  v = -0.5 / a[,"e2"]
  r = cbind(m = a[,"e1"] * v, v = v)
  # to avoid NaN if variance is infinite, set mean to something arbitrary
  r[ r[,"v"] == Inf, "m" ] = 0
  r
}

# Converts from canonical parameterization to mean and precision.
canonical.to.mean.and.precision = function(a) {
  p = -0.5 * a[,"e2"]
  cbind(m = a[,"e1"] / p, p = p)
}

# Log-evidence contribution for a normal variable.
# again, the question arises about whether "n" is needed...
# Args:
#   a - a list of normal distributions, as matrices with columns
#     "e1" and "e2"
# Returns: the log evidence contribution for those variables.
# FIXME not yet working
normal.log.evidence = function(a) {
  n = length(a)    # number of things being multiplied

  # convert to mean and precision
  a1 = sapply(a, canonical.to.mean.and.precision)

  # the product of these comes from summing parameters as usual
  s = 0 * a[[1]]
  for(i in 1:n)
    s = s + a[[i]]
  s1 = canonical.to.mean.and.precision(s)

  # add up precision portion
  sum1 = - sum( log(s1[,"p"]) )   # add "na.rm = TRUE" ?
#  for(i in 1:n)
#    sum1 = sum1 + sum(log(a1[[i]][,"p"]

#  sqrt( (p1 * p2 / p) / (2*pi) ) * exp(-0.5 * (p1 * m1^2 + p2 * m2^2 - p * m^2))

#  sum2 = - sum(
#  0.5 * (sum1 + FIXME) - 0.5 * sum2
}

# An equality factor for a normal variable.
# ??? is this right? I assume this is the integral
# along the line y = x.
# ETA: sort of. This should actually moment-match the
# product of the distributions.
normal.eq.factor = function(m) {
  a1 = canonical.to.mean.and.variance(m$a)
  b1 = canonical.to.mean.and.variance(m$b)



  list(m = list(a=m$b, b=m$a), log.evidence = FIXME)
}

# A normal prior factor.
normal.prior.factor = function(b) list(
  # this should force the outgoing message to be b
  # ??? is this right?
  update = function(a) list(x = a$x + b),
  log.evidence = function(a) {
    0   # FIXME
  }
)


