# Factor which forces a variable to be positive.

# Moment-matches a normal, truncated at x >= 0,
# with mean and variance m and v, respectively.
# following Phoebus Dhrymes'
# "Moments of Truncated (Normal) Distributions".
positive.moment.match = function(m, v) {
  m1 = -as.vector(m)
# v[v<0] = 1e10   # XXX hack
  s = sqrt(v)
  z = -m1 / s
  a = dnorm(z) / pnorm(z)

  # hack to deal with when z is very negative
#  r = cbind(m = ifelse(z < -30, 0, - (m1 - s * a)),
#    v = ifelse(z < -30, 0, v * (1 - z*a - a^2)))
  r = cbind(m = - (m1 - s * a), v = v * (1 - z*a - a^2))
  r
}


# Computes the message and log-evidence.
positive.factor = function(a) {
  b = canonical.to.mean.and.variance(a$x)
  r = positive.moment.match(b[,"m"], b[,"v"])
  list(
    m = list(x = mean.and.variance.to.canonical(r)),
    log.evidence = sum(pnorm(0, -b[,"mu"], sd, log.p=TRUE)))
}

