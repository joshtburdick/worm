# Normally-distributed variables, and related utilities.

# Creates a multivariate normal vector variable.
norm.var = function(n) {
  b = matrix(0, nrow=n, ncol=3)
  # XXX I didn't know what to call these, but named them this, since
  # Wikipedia calls the canonical parameters "eta"
  # ??? is "n" needed?
  colnames(b) = c("n", "e1", "e2")
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

