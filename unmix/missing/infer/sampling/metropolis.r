# Metropolis sampling.

library(mnormt)

# Metropolis update rule.
metropolis = function(x, x1, ll) {
  a = ll(x)
  b = ll(x1)
  (b >= a) || rbinom(1, 1, exp(b-a))
}

# Convenience function which actually does
# the update.
metropolis.update = function(x, x1, ll) {
  if (metropolis(x, x1, ll))
    x1
  else
    x
}

# A quick test: samples from a multivariate normal distribution.
metropolis.test.1 = function(mu, si, num.samples) {
  n = length(mu)
  x = rnorm(n)

  f = function(x)
    dmnorm(x, mean=mu, varcov=si, log=TRUE)

  X = matrix(nrow=num.samples, ncol=n)
  for(i in 1:num.samples) {
    X[i,] = x
    x1 = as.vector(rmnorm(mean=x, varcov=diag(n)))
    x = metropolis.update(x, x1, f)
  }

  X
}

metropolis.test.1a = function() {
  mu = c(3,4)
  si = matrix(c(1,0.8,0.8,1), nrow=2)

  X = metropolis.test.1(mu, si, 10000)
  plot(X[,1], X[,2], pch=20, col="#00000040")
}


