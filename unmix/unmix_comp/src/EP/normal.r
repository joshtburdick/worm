# Normally-distributed variables, and related utilities.

# Converts from mean and variance to canonical parameters.
mean.and.variance.to.canonical = function(a) {
  cbind(e1 = a[,"m"] / a[,"v"], e2 = -0.5 / a[,"v"])
}

# Converts from canonical parameters to mean and variance.
canonical.to.mean.and.variance = function(a) {
  v = -0.5 / a[,"e2"]
  r = cbind(m = a[,"e1"] * v, v = v)
  # to avoid NaN if variance is infinite, set mean to something arbitrary
  r[ r[,"v"] == Inf, "m" ] = 0
  r
}

