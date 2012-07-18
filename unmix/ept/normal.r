# Normally-distributed variables, and related utilities.



# Creates a multivariate normal vector variable.
norm.var = function(n) {
  b = matrix(0, nrow=3, ncol=n)
  # XXX I didn't know what to call these, but named them this, since
  # Wikipedia calls the canonical parameters "eta"
  rownames(b) = c("n", "e1", "e2")
  list(b = b, observed=FALSE)
}

# FIXME: these are not yet working

# Converts from canonical parameterization to mean and s.d.
canonical.to.mean.and.sd = function(a) {
  p = -2 * (a[,"e2"] / a[,"n"])
  cbind(n = rep(1, length(p)),
    mu = a[,"e1"] / a[,"n"],
    sd = sqrt(1/p))
}

# Converts from mean and s.d. to canonical parameters.
mean.and.sd.to.canonical = function(a) {
  p = 1 / (a[,"sd"]^2)
  cbind(n = rep(1, length(p)),
    e1 = a[,"n"] * (a[,"m"] / (a[,"sd"]^2)),
    e2 = -0.5 * a[,"n"] * p)
}


# Factor constraining some variables to be linear transforms of each other.








