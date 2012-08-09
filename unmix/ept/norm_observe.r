
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


# Adds an observation to a multivariate normal.
mvnorm.constrain = function(A, b, b.var) list(

  # update conditions on this observation
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



