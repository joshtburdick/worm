# A "scaled Dirichlet" ("SD") distribution, and various
# utilities. These aren't exact, but are fairly simple.

source("git/unmix/ept/gamma.r")
source("git/unmix/ept/moment.r")

# Converts from SD parameters to marginal mean and variance.
sd.s2mv = function(a) {
  mv = gamma.s2mv(a)
  moments = mv2moment(mv)
  s = apply(mv, 1, sum)
  s.moments = c(mv2moment(cbind(s)))
  moment2mv( moments / s.moments )
}

# Converts from marginal mean and variance to SD parameters.
sd.mv2s = function(b) {
  # FIXME: make sure these add up to 1?

  m = mv2moment(b)

  a = sum( m["x2",] )
  b = sum( m["x1",] ^ 2 )

  # the first moments stay the same
  m["x2",] = m["x2",] * ( (1-b) / (1-a) )

  gamma.mv2s( moment2mv(m) )
}

# some simple tests. For concreteness, I'm scaling
# the a's to sum to one (even though it doesn't affect
# the result of sd.s2mv()).
sd.test1 = function() {
  a1 = rbind(a=c(1,1,1), b=c(1,1,1))
  a1["b",] = a1["b",] * sum(gamma.s2mv(a1)["m",])
  b1 = sd.s2mv(a1)

  a2 = rbind(a=c(1,2,3), b=c(3,2,1))
  a2["b",] = a2["b",] * sum(gamma.s2mv(a2)["m",])
  b2 = sd.s2mv(a2)
  list(a1, b1, a2, b2)
}

# More simple tests.
sd.test2 = function(m, n) {
  for(i in 1:n) {
    a = rbind(a = rgamma(m,1,1), b = rgamma(m,1,1))
    b = sd.s2mv(a)
    b1 = sd.s2mv( sd.mv2s(b) )
    cat("diff =", range(b1 - b), "\n")
  }
}

sd.sampling.test.1 = function(n) {

  # for testing
  source("git/unmix/ept/practice/gamma_conditional_6.r")


  # standard parameters of some random gamma distributions
  s = rbind(a = sample(c(1:10), n, replace=TRUE),
    b = sample(c(1:5), n, replace=TRUE))
#    b = rep(1, n))
  print(s)
  mv.sd = sd.s2mv(s)
  print(mv.sd)
  mv.sampling = gamma.n2mv(gamma.cond.sampling(gamma.s2n(s), t(rep(1,n)), 1))
  print(mv.sampling)
  cat("\n")
}
# these are pretty convincingly different

