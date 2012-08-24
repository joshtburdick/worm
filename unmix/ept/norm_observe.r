
# Adds an observation to a normal variable.
# Args:
#   m - mean of observation
#   sd - standard deviation (if this is a "soft" observation)
norm.observe = function(m, sd = 0) list(

  update = function(   ) {

  },

  log.evidence = function(   ) {

  }
)

# Distribution of x ~ N(m,v) | Ax ~ N(b,b.var).
# FIXME: this doesn't deal well with when any part of v is infinite
# (although that seems to be well-defined.)
lin.constraint.1 = function(m, v, A, b, b.var) {
  M = A %*% (v * t(A))
  diag(M) = diag(M) + b.var

  # using Cholesky would be slightly faster, but maybe less stable
  M.qr = qr(M)

  r = cbind(m = m - v * as.vector(t(A) %*% solve(M.qr, A %*% m - b)),
    v = v - v * apply(A * solve(M.qr, A), 2, sum) * v)
  r
}

# Distribution of x ~ N(m,v) | Ax ~ N(b,b.var).
# Like lin.constraint.1(), but with "v" possibly infinite.
lin.constraint.2 = function(m, v, A, b, b.var) {

  inf.v = is.infinite(v)

  # if there are strictly more entries with v infinite than there are
  # constraints, then it seems like the posterior covariance is singular.
  if (sum(inf.v) > nrow(A)) {
    # entries with finite variance are simply copied through.
    # entries with infinite variance are set to a large arbitrary number
    m[inf.v] = 0
    v[inf.v] = 1e10
  }
  else {
    # if there are at least as many constraints as there are infinite
    # entries of v, then the posterior should be well-defined.




  }

  # get finite entries of v, and those elements of A
  v1 = v[ is.finite(v) ]
  A1 = A[ , is.finite(v) ]
  m1 = m[ is.finite(v) ]

  M = A1 %*% (v1 * t(A1))
  diag(M) = diag(M) + b.var

  # using Cholesky would be slightly faster, but maybe less stable
  M.qr = qr(M)

  m1 - v1 * as.vector(t(A1) %*% solve(M.qr, A1 %*% m1 - b))

#  r = cbind(m = m - v * as.vector(t(A) %*% solve(M.qr, A %*% m - b)),
#    v = v - v * apply(A * solve(M.qr, A), 2, sum) * v)
#  r
}


# Adds an observation to a multivariate normal.
mvnorm.constrain = function(A, b, b.var) list(

  # update conditional on this observation
  update = function(a) {
    mv = canonical.to.mean.and.variance(a$x)
    r = lin.constraint.1(mv[,"m"], mv[,"v"], A, b, b.var)
    mean.and.variance.to.canonical(r)
  },

  # this is basically the likelihood for this observation
  log.evidence = function(a) {
    0    # FIXME
  }
)

mvnorm.observe = function(A, b, b.var) function(a) {




}

