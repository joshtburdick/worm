# Converting to and from moments of the inverse-gamma function,
# in order to compute marginals of this thing.

source("git/unmix/ept/gamma.r")
source("git/unmix/ept/moment.r")

# Computes mean and variance of an inverse gamma.
# (Note that the mean is only defined when a > 1,
# and the mode is only defined when a > 2.)
# Args:
#   a - parameters of the inverse gamma (specifically,
#     if X ~ Gamma(a, b), where a and b are shape and
#     rate respectively, then this has distribution 1/X)
# Returns: mean and variance of that distribution
ig.s2mv = function(a) {
  rbind(m = a["b",] / (a["a",] - 1),
    v = (a["b",] ^ 2) / (((a["a",] - 1)^2) * (a["a",] - 2)))
}

# Computes moments of X/(X+Y), where X and Y have
# gamma distributions.
# Args:
#   x, y - standard parameters of the relevant gamma
#     distributions
# Returns: moments of X/(X+Y)
gamma.ratio.moments = function(x, y) {

  # first, moments of X/Y
  x.m = mv2moment( gamma.s2mv(x) )
  y1.m = mv2moment( ig.s2mv(y) )
  a = moment2mv(x.m * y1.m)
  print("x/y")
  print(a) 
  a["m",1] = a["m",1] + 1
  print("1 + x/y")
  print(a)

  print("parameters =")
  a.p = gamma.mv2s(a)
  print(a.p)

  # compute estimate of 1/(all that)
  a1 = ig.s2mv( a.p )
  a1
}

# Test of gamma.ratio.moments().
gamma.ratio.moments.test = function(xa, xb, ya, yb) {

  # first, test estimate
  x1 = rbind(a=xa, b=xb)
  y1 = rbind(a=ya, b=yb)
  print(gamma.ratio.moments(x1,y1))

  # then, test by sampling
  n = 1e7
  x = rgamma(n, shape=xa, rate=xb)
  y = rgamma(n, shape=ya, rate=yb)
# browser()
  mean.and.var = function(a) {
    paste("mean =", round(mean(a),4), "  var =", round(var(a),4))
  }
  cat("sampling:\n")
  cat(paste0("x/y: ", mean.and.var(x/y), "\n"))
  cat(paste0("1+x/y: ", mean.and.var(1+x/y), "\n"))
  # XXX note that 1+(x/y) is always greater than 1; this implies that
  # a gamma distribution will be a bad match for it

  cat(paste0("1/(1+x/y): ", mean.and.var(1/(1+x/y)), "\n"))
}

# Simple tests of the inverse gamma.
ig.test1 = function() {
  a1 = rbind(a=c(3,4,5), b=c(2,5,1))
  print(ig.s2mv(a1))
  x1 = 1 / rgamma(1e6, shape=3, rate=2)
  print(mean(x1))
  print(var(x1))
  x1 = 1 / rgamma(1e6, shape=4, rate=5)
  print(mean(x1))
  print(var(x1))
  x1 = 1 / rgamma(1e6, shape=5, rate=1)
  print(mean(x1))
  print(var(x1))
}

gamma.ratio.moments.test(3,1,4,1)




