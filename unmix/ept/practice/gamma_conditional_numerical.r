# Moments of gamma-distributed numbers summing to 1, using
# quadrature.

source("git/unmix/ept/gamma.r")
source("git/unmix/ept/practice/gamma_conditional_6.r")

# Moments of two gamma-distributed vars,
# conditional on them summing to 1.
# Args:
#   g, h - standard parameters of two gamma distributions
#   x - where to evaluate the density
# Returns: mean and variance, based on those points
gamma.cond.sum.numerical.int.1 = function(g, h) {
  knots = c(1:10000) / 10000
  g = as.vector(g)
  h = as.vector(h)
  f = dgamma(knots, shape = g[1], rate = g[2], log = TRUE) +
    dgamma(1 - knots, shape = h[1], rate = h[2], log = TRUE)
  f = f - max(f)
  w = exp(f)
  w = w / sum(w)
  w[ is.na(w) ] = 0
  w[ w < 1e-300 ] = 0
  x1 = sum(w * knots)
  x2 = sum(w * (knots^2))
  c(m = x1, v = x2 - x1^2)
}

# like gamma.cond.sum.numerical.int.1, but uses built-in integration.
# FIXME not working
gamma.cond.sum.numerical.int.2 = function(g, h) {
  f = function(x) dgamma(x, shape = g[1], rate = g[2]) *
    dgamma(1 - x, shape = h[1], rate = h[2])
  x0 = integrate(function(x) f(x), 0, 1)$value
  x1 = integrate(function(x) f(x) * x / x0, 0, 1)$value
  x2 = integrate(function(x) f(x) * (x^2) / x0, 0, 1)$value

  c(m = x1, v = x2 - x1^2)
}

# Moments of gamma-distributed vars, conditional on
# them summing to 1.
gamma.cond.sum.numerical.1 = function(x) {
  g = gamma.n2s(x)
  mv = gamma.n2mv(x)
  s = apply(mv, 1, sum)
  mv.other = s - mv

  r = 0 * mv
  g.other = gamma.mv2s(mv.other)

  for(j in 1:ncol(r)) {
    r[,j] = gamma.cond.sum.numerical.int.1(g[,j], g.other[,j])
  }

  r
}

# Moments of gamma-distributed vars, conditional on
# a linear constraint.
gamma.cond.sum.numerical = function(x, a, b) {

 # convert to moments
  m = gamma.n2mv(x)

  # scale this
  s = b / a
  m["m",] = m["m",] / s
  m["v",] = m["v",] / (s*s)

  # approximate marginals, if those sum to 1
  r = gamma.cond.sum.numerical.1(gamma.mv2n(m))

  # undo scaling
  r["m",] = r["m",] * s
  r["v",] = r["v",] * (s*s)

  # convert back to natural parameters
  gamma.mv2n(r)
}

# some simple tests of the unscaled version
if (FALSE) {
  # Test of this, comparing to sampling.
  # Args: x - standard params of gamma distributions
  gamma.cond.sum.numerical.t1 = function(x) {
    x.n = gamma.s2n(x)

    cat("estimated:\n")
    print(gamma.cond.sum.numerical.1(x.n))

    cat("from sampling:\n")
    print(gamma.n2mv(gamma.cond.sampling(gamma.s2n(x),
      t(rep(1, ncol(x))), 1)))

    cat("\n")
  }

  gamma.cond.sum.numerical.t1(rbind(a=c(1,1,1), b=c(1,1,1)))
  gamma.cond.sum.numerical.t1(rbind(a=c(1,2,3), b=c(1,1,1)))
  gamma.cond.sum.numerical.t1(rbind(a=c(1,1,1), b=c(3,1,1)))
  gamma.cond.sum.numerical.t1(rbind(a=c(2,5,9), b=c(4,1,1)))
}
# these are not exact, but are the closest I've seen
# (often variance from sampling is smaller, which is consistent
# with sampling not mixing completely)

# ... and tests of the scaled version
if (FALSE) {

  gamma.cond.sum.numerical.t2 = function(x, a, b) {
    x.n = gamma.s2n(x)

    cat("estimated:\n")
    print(gamma.n2mv(gamma.cond.sum.numerical(x.n, a, b)))

    cat("from sampling:\n")
    print(gamma.n2mv(gamma.cond.sampling(x.n, a, b)))

    cat("\n")
  }
  f = gamma.cond.sum.numerical.t2
  x = rbind(a=c(1,1,1), b=c(1,1,1))
#  f(x, t(c(1,1,1)), 1)
#  f(x, t(c(1,2,3)), 1)
  x = rbind(a=c(2,5,9), b=c(4,1,1))
  f(x, t(c(1,1,1)), 4)
  f(x, t(c(1,2,3)), 7)

  x = rbind(a=sample(1:10, 10, replace=TRUE), b=sample(1:10, 10, replace=TRUE))
  f(x, t(sample(1:10,10, replace=TRUE)), sample(1:5, 1))
}



