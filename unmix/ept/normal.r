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
  cbind(e1 = a[,"m"] / a[,"v"], e2 = -0.5 / a[,"v"])
}

# Converts from canonical parameterization to mean and variance.
canonical.to.mean.and.variance = function(a) {
  v = -0.5 / a[,"e2"]
  r = cbind(m = a[,"e1"] * v, v = v)
  # to avoid NaN if variance is infinite, set mean to something arbitrary
  r[ r[,"v"] == Inf, "m" ] = 0
  r
}

# Log-evidence contribution for a normal variable.
# again, the question about whether "n" is needed...
# Args:
#   a - a list of normal distributions, as matrices with columns
#     "n", "e1", and "e2"
# Returns: the log evidence contribution for those variables.
normal.log.evidence = function(a) {


  # does this just depend on the variance? No.
}

# An equality factor for a normal variable.
# ??? is this right? I assume this is the integral
# along the line y = x.
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


