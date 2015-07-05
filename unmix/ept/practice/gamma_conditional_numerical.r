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
gamma.cond.sum.numerical.1 = function(g, h) {
  knots = c(1:100000) / 100000
  g = as.vector(g)
  h = as.vector(h)
  f = dgamma(knots, shape = g[1], rate = g[2], log = TRUE) +
    dgamma(1 - knots, shape = h[1], rate = h[2], log = TRUE)
  f = f - max(f)
  w = exp(f)
  w = w / sum(w)
  x1 = sum(w * knots)
  x2 = sum(w * (knots^2))
  c(m = x1, v = x2 - x1^2)
}

# Moments of gamma-distributed vars, conditional on
# them summing to 1.
gamma.cond.sum.numerical = function(x) {
  g = gamma.n2s(x)
  mv = gamma.n2mv(x)
  s = apply(mv, 1, sum)
  mv.other = s - mv

  r = 0 * mv
  g.other = gamma.mv2s(mv.other)

  for(j in 1:ncol(r)) {
    r[,j] = gamma.cond.sum.numerical.1(g[,j], g.other[,j])
  }

  r
}

# some simple tests
if (FALSE) {
  # Test of this, comparing to sampling.
  # Args: x - standard params of gamma distributions
  gamma.cond.sum.numerical.t1 = function(x) {
    x.n = gamma.s2n(x)

    cat("estimated:\n")
    print(gamma.cond.sum.numerical(x.n))

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


