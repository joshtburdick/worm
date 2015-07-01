# Attempt at computing these marginals with the exact numbers,
# using likelihood weighting, and comparing these with
# approximations thereof.

source("git/unmix/ept/gamma.r")

library(MASS)

wd = getwd()
setwd("~/gcb/git/unmix/unmix_comp/src/")
source("sampling/cdaCpp.r")
setwd(wd)

xlim=c(0,5)


# Samples from Ax=b, x>=0.
sample.lin.pos = function(A, b) {
  n = ncol(A)

  x0 = ldei(E=A, F=b, G=diag(n), H=rep(0, n))$X

  sample.cda(A, b, x0, 1e6, n)   # was 1e6 
}

# Computes stats of weighted rows of matrices.
weighted.stats = function(x, w) {
  x0 = sum(w)
  x1 = apply(w * x, 2, sum) / x0
  x2 = apply(w * x^2, 2, sum) / x0

  rbind(m = x1,
    v = x2 - x1^2)
}

# Weights some samples.
# Args:
#   x - the samples
#   p - the gamma parameters (as natural parameters)
# Returns: list containing
#   w - the weights
#   p1 - the gamma parameters (conditional on x)
gamma.likelihood.weight = function(x, p) {
  p = gamma.n2s(p)
  stopifnot(ncol(x) == ncol(p))

  # add up (log-transformed) weights
  lw = rep(0, nrow(x))
  for(j in 1:ncol(x)) {
    lw = lw + dgamma(x[,j], shape=p["a",j], rate=p["b",j], log=TRUE)
  }

  # scale, and convert from log
  w = exp( lw - max(lw) )

  list(w = w, p1 = gamma.mv2n(weighted.stats(x, w)))
}

# Computes marginals of a gamma distribution by sampling.
# Args:
#   g - parameters of gamma distributions
#   A, b - these give the linear constraints
# Returns: list containing
#   x - the samples
#   g1 - parameters of marginals, conditional on Ax=b
#     (as natural parameters.)
gamma.cond.sampling = function(g, A, b) {
  x = sample.lin.pos(A, b)

  gamma.likelihood.weight(x, g)$p1
}

# Alternative (hopefully faster) version of that, in the
# case that A just has one row.
# Args:
#   g - parameters of gamma distributions
#   A, b - the linear constraint
# Returns: list (as from gamma.cond.sampling)
# Deprecated; not sure it works.
# cache of Gamma(1,1) samples.
# g.cache = matrix(rgamma(1e7,1,1), nrow=10)
gamma.cond.sampling.2 = function(g, A, b) {
  x = t( g.cache[ c(1:ncol(g)) , ] * as.vector(1/A) )
  x = b * x / apply(x, 1, sum)

  gamma.likelihood.weight(x, g)$p1
}

# x1 = sample.lin.pos(t(c(1,1,1,1,1)), 1)

# comparing marginals for different orders of conditioning things
x2 = gamma.s2n(rbind(a=rep(1,10), b=rep(1,10)))
# A2 = matrix(c(1,1,0.1, 1,0.1,1), nrow=2, byrow=TRUE)
# b2 = c(2,3)
A2 = matrix(rgamma(20, shape=1, rate=1), nrow=2, byrow=TRUE)
y2 = rgamma(10, shape=0.5, rate=1)
b2 = as.vector(A2 %*% y2)

# compares conditioning using one constraint Ax = b, with
# an approximation from adding each constraint successively
test1 = function() {
# not using this
# gamma.cond.sampling = gamma.cond.sampling.2
print(A2)
print(y2)
print(b2)

x2a = gamma.cond.sampling(x2, A2, b2)
print(gamma.n2mv(x2a))

x2b1 = gamma.cond.sampling(x2, t(A2[1,]), b2[1])
x2b2 = gamma.cond.sampling(x2b1, t(A2[2,]), b2[2])
print(gamma.n2mv(x2b2))

x2c1 = gamma.cond.sampling(x2, t(A2[2,]), b2[2])
x2c2 = gamma.cond.sampling(x2c1, t(A2[1,]), b2[1])
print(gamma.n2mv(x2c2))

par(mfrow=c(2,3))
plot(gamma.n2mv(x2a)[1,], gamma.n2mv(x2b2)[1,])
abline(0,1)
plot(gamma.n2mv(x2a)[1,], gamma.n2mv(x2c2)[1,])
abline(0,1)
plot(gamma.n2mv(x2b2)[1,], gamma.n2mv(x2c2)[1,])
abline(0,1)

plot(gamma.n2mv(x2a)[2,], gamma.n2mv(x2b2)[2,])
abline(0,1)
plot(gamma.n2mv(x2a)[2,], gamma.n2mv(x2c2)[2,])
abline(0,1)
plot(gamma.n2mv(x2b2)[2,], gamma.n2mv(x2c2)[2,])
abline(0,1)

}

# these are in fact different, even with the "exact" computation

# XXX but this seems incorrect, as e.g. this doesn't have a mean adding up to 1:
# oh, wait, now it does.
# gamma.n2mv(gamma.cond.sampling(gamma.s2n(rbind(a=rep(1,3), b=rep(1,3))), t(c(1,1,1)), 1))

# test1()


