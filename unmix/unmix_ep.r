# Unmixes using EP.

source("git/unmix/ept/normal.r")
# source("git/unmix/ept/matrix_inv_lemma.r")   XXX um, don't think I need this

# Moment-matches a normal, truncated at x >= 0,
# with mean and variance m and v, respectively.
# following Phoebus Dhrymes'
# "Moments of Truncated (Normal) Distributions".
norm.moment.match = function(m, v) {
  m1 = -as.vector(m)
v[v<0] = 1e8   # XXX hack
  s = sqrt(v)
  z = -m1 / s
  a = dnorm(z) / pnorm(z)

  # hack to deal with when z is very negative
  r = cbind(m = ifelse(z < -30, 0, - (m1 - s * a)),
    v = ifelse(z < -30, 0, v * (1 - z*a - a^2)))
  r
}

# Messages constraining a variable to be positive.
positive.factor = function(x) {
  mv = canonical.to.mean.and.variance(x)
  mm = norm.moment.match(mv[,"m"], m.and.v[,"v"])
  mean.and.variance.to.canonical(mm)
}

# Distribution of x ~ N(m,v) | Ax ~ N(b,b.var).
lin.constraint = function(m, v, A, b, b.var) {

  # first, compute (A V A^T + diag(v)) ^ -1, which I'll call B
  # FIXME: the first Diagonal() call here can probably be simplified
  M = as.matrix( A %*% Diagonal(x=V) %*% t(A) + Diagonal(x = b.var) )
  B.chol = chol2inv(chol(B))

  cbind(m = m - v * t(A) %*% solve(B.chol, A %*% m - b),
    v = v - v * apply(t(A) * solve(B.chol, A), 1, sum) * v)
}

# Message constraining "Ax ~ N(b, b.var)".
lin.constraint.factor = function(A, b, b.var) function(x) {
  mv = canonical.to.mean.and.variance(x)
  r = lin.constraint(mv[,"m"], mv[,"v"], A, b, b.var)
  mean.and.variance.to.canonical(r)
}

# Does unmixing when the constraints are
#   Ax ~ Normal(-,-), x >= 0.
# (Assumes a flat prior.)
unmix.ep = function(A, b, b.var) {
  n = length(b)

  # the term approximations (initially arbtrary, with nonzero variance)
  terms = mean.and.variance.to.canonical(cbind(m=rep(0,n), v=rep(1,n)))

  # the posterior
  q = terms

  for(iter in 1:5) {



    # add in Ax ~ N(-,-) constraint



  }

  mv = canonical.to.mean.and.variance(q)
  list(m = mv[,"m"], v = mv[,"v"])
}

