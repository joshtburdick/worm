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

# A factor which forces a variable to be positive.
positive.factor = list(

  # moment-matches the positive portion
  update = function(a) {
    b = canonical.to.mean.and.variance(a)
    r = positive.moment.match(b[,"m"], b[,"v"])
    mean.and.variance.to.canonical(r)
  },

  # log of the proportion of the distribution that's positive
  log.evidence = function(a) {
    b = canonical.to.mean.and.variance(x)
    sum(pnorm(0, -b[,"mu"], sd, log.p=TRUE))
  }
)

