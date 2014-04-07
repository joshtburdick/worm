# Approximates a constrained region using EP --
# version designed to deal with a full covariance matrix.

library(Matrix)
library(corpcor)

source("git/utils.r")
source("git/unmix/ept/gamma.r")
source("git/unmix/ept/normal.r")

# Moment-matches a normal, truncated at x >= 0,
# with mean and variance m and v, respectively.
# following Phoebus Dhrymes'
# "Moments of Truncated (Normal) Distributions".
# This version doesn't check for anything being
# negative.
positive.moment.match = function(m, v) {
  m1 = - as.vector(m)
  s = sqrt(v)
  z = -m1 / s
  a = dnorm(z) / pnorm(z)

  r = cbind(m = - (m1 - s * a), v = v * (1 - z*a - a^2))
}

# Messages constraining a variable to be positive
# (but using canonical parameterization.)
positive.moment.match.canonical = function(x) {
  mv = canonical.to.mean.and.variance(x)
  mm = positive.moment.match(mv[,"m"], mv[,"v"])
  mean.and.variance.to.canonical(mm)
}

# Adds in a constraint of the form Ax ~ b, and returns the
# full covariance matrix. This expects, and returns, a
# full covariance matrix.
mvnorm.naive = function(m, V, A, b, b.var) {

  # first, compute (A V A^T + v I) ^ -1, which I'll call B
  M = as.matrix( A %*% V %*% t(A) + Diagonal(x = b.var) )

  B = pseudoinverse(M)

  list(m = m - as.vector(V %*% t(A) %*% B %*% (A %*% m - b)),
    V = as.matrix( V - V %*% t(A) %*% B %*% A %*% V) )
}

# Distribution of x ~ N(m,v) | Ax ~ N(b,b.var).
lin.constraint.1 = function(m, v, A, b, b.var) {

  # first, compute (A V A^T + diag(v)) ^ -1, which I'll call B
#  M = as.matrix( A %*% (v * t(A)) + Diagonal(x = b.var) )
  M = A %*% (v * t(A))
  diag(M) = diag(M) + b.var

  M.pi = pseudoinverse(M)

  # trying to use the pseudoinverse here
  r = cbind(m = m - v * as.vector(t(A) %*% (M.pi %*% (A %*% m - b))),
    v = v - v * apply(A * (M.pi %*% A), 2, sum) * v)

  r
}

# Distribution of x ~ N(m,v) | Ax ~ N(b,b.var),
# allowing A to be sparse.
lin.constraint.sparse = function(m, v, A, b, b.var) {

  # first, compute (A V A^T + diag(v)) ^ -1, which I'll call B
#  M = as.matrix( A %*% (v * t(A)) + Diagonal(x = b.var) )
  M = A %*% (v * t(A))
  diag(M) = diag(M) + b.var

  M.chol = chol(M)

  r = cbind(
    m = m - v * as.vector(t(A) %*% chol.solve(M.chol, A %*% m - b)),
    v = v - v * colSums(A * chol.solve(M.chol, A)) * v)
  r
}

# Message constraining "Ax ~ N(b, b.var)".
lin.constraint.factor = function(A, b, b.var) function(x) {
  mv = canonical.to.mean.and.variance(x)
#  r = mvnorm.2.diag(mv[,"m"], mv[,"v"], A, b, b.var)
  r = lin.constraint.sparse(mv[,"m"], mv[,"v"], A, b, b.var)
# print("r =")
# print(r)
# print(colnames(r))
  mean.and.variance.to.canonical(r)
}

# Approximates a region constrained such that x >= 0,
# allowing a full covariance matrix.
# Args:
#   A, b, b.var - these give the constraint that
#     Ax ~ Normal(b, b.var), x >= 0
#     (For exact constraints, b.var can be 0.)
#   prior.mean - the prior mean (zero if omitted)
#   prior.var - the prior covariance (flat if omitted)
# Returns: list with components
#   m, v - the mean and variance of the posterior
#   t - the prior times the terms (without the linear constraint)
#   update.stats - matrix with mean and variance of update sizes
approx.region.full = function(A, b, b.var,
  prior.mean=NULL, prior.var=NULL, converge.tolerance = 1e-9, max.iters=100) {
  n = ncol(A)

#  debug.dir = "~/tmp/approx.region.debug"
  debug.dir = NULL

  # the term approximations (initially flat)
  terms = mean.and.variance.to.canonical(cbind(m=rep(0,n), v=rep(Inf,n)))

  # the (possibly soft) linear constraint
  lin.constraint = lin.constraint.factor(A, b, b.var)

  # the posterior; note that we only track the diagonal elements of this
  # ??? is this correct?
  q = mean.and.variance.to.canonical(cbind(m=rep(0,n), v=rep(1,n)))
#  q = lin.constraint( prior + terms )
#  q = prior
  # convergence statistics
  update.stats = NULL

  for(iter in 1:max.iters) {
    q.old = q
    terms.old = terms

    terms.1 = q - terms

    mm = positive.moment.match.canonical(terms.1)

    # update terms (currently without any damping)
    terms = (mm - q) + terms

    if (!is.null(debug.dir)) {
      system(paste("mkdir -p", debug.dir))
      ep.trace = list(terms=terms, terms.1=terms.1, q=q, mm=mm)
      save(ep.trace, file=paste(debug.dir, "/", iter, ".Rdata", sep=""))
    }

    # add in Ax ~ N(-,-) constraint
    q = lin.constraint( prior + terms )

    # ??? show change in mean and variance separately?
    diff = apply(abs(canonical.to.mean.and.variance(q) - canonical.to.mean.and.variance(q.old)), 2, max)
 write.status(paste(signif(diff, 2)))
    update.stats = rbind(update.stats, diff)

    # possibly stop early
    if (max(diff, na.rm=TRUE) <= converge.tolerance)
      break
  }
cat(" \n")

  mv = canonical.to.mean.and.variance(q)
  list(m = mv[,"m"], v = mv[,"v"],
    t = canonical.to.mean.and.variance( prior + terms ),
    update.stats = update.stats)
}

