# Attempt at implementing Minka's clutter EP algorithm.
# For now, attempting a very literal (but slow and awkward)
# transliteration of the math in the paper, mostly from
# section 3.1 .

library(mnormt)

# Creates the initial model (as in section 3.1).
# Note that all "scale factors" s are log-transformed.
init.model = function(n, d, w) {

  # prior variance
   v.0 = 100

  list(

    # number of examples, dimensionality, and weight of the
    # "signal" peak with unknown mean (as opposed to background)
    n = n, d = d, w = w,

    # current estimate of the "signal" peak (starts at the prior)
    m.x = rep(0, d), v.x = v.0,

    # the prior term
    m.0 = rep(0, d), v.0 = v.0, s.0 = (-d/2) * log(2*pi*v.0),

    # the data terms
    m = matrix(0, nrow=d, ncol=n), v = rep(1e4, n), s = rep(0, n))
}

# Density of a multivariate normal with spherical covariance.
# Deprecated (using dmnorm instead.)
mvnorm.sphere = function(x, m, v) {
  if (v <= 0) {
    cat("in mvnorm, v = ", v, "\n")
    v = 10^8     # XXX
  }
  exp( (-1/(2*v)) * sum((x-m)^2) )
}

# Does one step of EP.
# ??? am I doing these in the right order?
ep.update = function(m, y) {

  for(i in 1:m$n) {

    # 3.(a): "old" posteriors for data terms (the "...\i" term)
    old.v = 1 / ( 1/m$v.x - 1/m$v[i] )
    old.m = m$m.x + (old.v / m$v[i]) * t( m$m.x - m$m[,i] )
# print(c(i, m$v.x, old.v, old.m))

    # 3.(b) (as in ADF): update estimates of m.x and v.x
    r = 1 - (1/z) * m$w * dmnorm(y[,i],
      rep(m$m.0, m$d), sqrt(m$v.0) * diag(m$d))
#    l.off = m$w * mvnorm.sphere(y[,i], rep(0, m$d), 10)
#    l.on = (1 - m$w) * mvnorm.sphere(y[,i], old.m, old.v + 1)
#    z = l.off + l.on
#    r = l.on / z
# print(c(i,r))
    m$m.x = old.m + old.v * r * (y[,i] - old.m) / (old.v + 1)
    m$v.x = old.v - r * (old.v^2) / (old.v+1) +
      r * (1-r) * (old.v^2) * sum((y[,i] - old.m)^2) / (m$d * ((old.v+1)^2))
    z_i = (1 - m$w) * dmnorm(y[,i], old.m, (old.v+1) * diag(m$d)) +
      m$w * dmnorm(y[,i], rep(0, m$d), m$v.0 * diag(m$d))

    # 3.(c): update this term approximation
    m$v[i] = 1 / (1 / m$v.x - 1 / old.v)
# print(c(i, m$v.x, vi.before.update, m$v[i]))
#    if (is.finite(m$v[i])) {
      m$m[,i] = old.m + ((m$v[i] + old.v) / old.v) * (m$m.x - old.m)
      m$s[i] = log(z_i) - (m$d/2) * log(2 * pi * m$v[i]) -
        dmnorm(m$m[,i], old.m, (m$v[i] + old.v[i]) * diag(m$d), log=TRUE)
#    }
#    else {    FIXME implement this
      # XXX if a data point has zero likelihood, ignore it
#      m$s[i] = 0
#    }
  }

  m
}

# Computes the (log of the) evidence (and normalizing constant -- step 4.)
ep.log.evidence = function(m, y) {

  B = t(m$m.x) %*% m$m.x / m$v.x
    - sum( t(m$m) %*% m$m / m$v )

  (m$d/2) * log(2 * pi * m$v.x) + (B / 2) + sum(m$s)
}


# Constructs a toy data set.
set.seed(19)
clutter.y = matrix(rnorm(200, sd=10), nrow=2)
signal.y = matrix(rnorm(200), nrow=2) + 5
y = cbind(clutter.y, signal.y)
# y = y[,sample(10)]

m = init.model(100, 2, 0.5)

m1 = ep.update(m, y)


