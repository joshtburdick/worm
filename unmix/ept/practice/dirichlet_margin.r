# Attempt at "observing a marginal of a Dirichlet."

library("gtools")

# Given a Dirichlet, updates it with a "hard" observation
# (or attempt at same.)
# Args:
#   alpha - Dirichlet parameters
#   m - number of elements in the first set
# Returns: Dirichlet, after observing that split
dirichlet.margin = function(alpha, m, p) {
  n = length(alpha)
  i1 = 1:m
  i2 = (m+1):n

  s = sum(alpha)
  a1 = sum(alpha[i1]) / s
  b1 = sum(alpha[i2]) / s

  alpha1 = alpha
  alpha1[i1] = alpha1[i1] * (p / a1)
  alpha1[i2] = alpha1[i2] * ((1-p) / b1)
  alpha1
}

# Given samples, estimates a Dirichlet with the same moments.
dirichlet.est.moments = function(x, w = NULL) {

  # default weighting
  if (is.null(w))
    w = rep(1, nrow(x))

  # (weighted) moments
  x1 = apply(w*x, 2, sum) / sum(w)
  x2 = apply(w*x*x, 2, sum) / sum(w)

  # different estimates of the concentration
  s = (x1 - x2) / (x2 - x1^2)

  # Minka mentions a suggestion from Ronning which may be
  # better, but I'll just take the geometric mean, basically.
  alpha = (exp(mean(log(s+1))) - 1) * x1

  list(x1 = x1, x2 = x2, s = s, alpha = alpha)
}

# Quick test of the above.
test1 = function(a = c(4,7,11,6)) {
  x = rdirichlet(1e6, a)
  m = dirichlet.est.moments(x)
  print(a)
  print(round(m$alpha, 2))
}

# Samples from a Dirichlet, then conditions on the first
# few elements having a particular Beta distribution.
sample.dirichlet.beta.1 = function(samples=1e6) {
  a = sample(1:10, 5)
  b = sample(1:10, 2)

  x = rdirichlet(samples, a)
  s = apply(x[,1:3], 1, sum)
  w = dbeta(s, b[1], b[2])
  x.param.estimate = dirichlet.est.moments(x, w)

  r = list(a=a, b=b, p = x.param.estimate$alpha, n = sum(w))
}

sample.many = function() {
  dirichlet.beta.samples = NULL
  for(i in 1:1000) {
cat(i, "")
    r1 = c(sample.dirichlet.beta.1(1e6), recursive=TRUE)
    dirichlet.beta.samples = rbind(dirichlet.beta.samples, r1)
    if (i %% 10 == 0)
      save(dirichlet.beta.samples,
        file="git/unmix/ept/practice/dirichlet.beta.samples.Rdata")
  }
}

