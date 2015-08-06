# Moments of gamma-distributed numbers summing to 1, using
# quadrature.

source("git/unmix/ept/gamma.r")
source("git/unmix/ept/moment.r")
source("git/unmix/ept/practice/gamma_conditional_6.r")

# Integration using the trapezoid rule.
trapezoid.rule = function(x, y) {
  n = length(x)
  0.5 * sum( (x[-1] - x[-n]) * (y[-1] + y[-n]) )
}

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
  w[ is.na(w) ] = 0
  w[ w < 1e-300 ] = 0
  w = w / sum(w)
  x1 = sum(w * knots)
  x2 = sum(w * (knots^2))
  c(m = x1, v = x2 - x1^2)
}

# like gamma.cond.sum.numerical.int.1, but uses built-in integration.
# FIXME not working so well
gamma.cond.sum.numerical.int.2 = function(g, h) {
  f = function(x) dgamma(x, shape = g[1], rate = g[2]) *
    dgamma(1 - x, shape = h[1], rate = h[2])
  x0 = integrate(function(x) f(x), 0, 1)$value
  x1 = integrate(function(x) f(x) * x / x0, 0, 1)$value
  x2 = integrate(function(x) f(x) * (x^2) / x0, 0, 1)$value

  c(m = x1, v = x2 - x1^2)
}

# Golden section search. Adapted from
# https://en.wikipedia.org/wiki/Golden_section_search,
# except that this maximizes instead.
golden.section.search = function(f, A, B, n=10) {
  gr = (sqrt(5)-1) / 2
  x = rep(NA, n)
  y = rep(NA, n)
  x[1] = A
  y[1] = f(A)
  x[2] = B
  y[2] = f(B)

  C = B - gr*(B-A)
  D = A + gr*(B-A)
  for(iter in 3:n) {
    fC = f(C)
    fD = f(D)
    if (fC > fD) {
      B = D
      D = C
      C = B - gr*(B-A)
      x[iter] = C
      y[iter] = f(C)
    }
    else {
      A = C
      C = D
      D = A + gr*(B-A)
      x[iter] = D
      y[iter] = f(D)
    }
  }

  list(x = x, y = y)
}


# Numerical way of computing this, using a grid of points
# from the golden section search.
# Args:
#   g, h - standard parameters of two gamma distributions
# Returns: natural parameters of the corresponding
#   moment-matched gamma distribution.
gamma.cond.sum.numerical.int.3 = function(g, h) {
  f = function(x)
    dgamma(x, shape = g[1], rate = g[2]) *
    dgamma(1 - x, shape = h[1], rate = h[2])
  r = golden.section.search(f, 0, 1, n = 100)
# browser()
#  y = exp( r$y - max(r$y) )
  w = r$y
  w[ is.na(w) ] = 0
  w[ w < 1e-300 ] = 0
  w = w / sum(w)

# XXX this way of finding moments doesn't seem to work
if (TRUE) {
  a = gamma.mv2n(moment2mv(
    rbind(x1 = sum(w * r$x), x2 = sum(w * (r$x ^ 2)) )))
  b = gamma.mv2n(moment2mv(
    rbind(x1 = mean(r$x), x2 = mean(r$x^2))))
  return(gamma.n2mv(a - b))
  # XXX not working
}

  i = order(r$x)
  x = r$x[i]
  y = r$y[i]
  n = length(x)
  s = trapezoid.rule(x, y)

# gamma.mv2n
  (moment2mv(
    rbind(x1 = trapezoid.rule(x, y * x) / s,
      x2 = trapezoid.rule(x, y * (x^2)) / s)))
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


foo = gamma.cond.sum.numerical.int.3(rbind(a=1,b=1), rbind(a=1,b=1))

