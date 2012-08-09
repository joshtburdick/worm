# Approximates a constrained region using EP.

library(Matrix)
library(corpcor)

source("git/unmix/ept/normal.r")
# source("git/unmix/ept/matrix_inv_lemma.r")   XXX um, don't think I need this

backspace = cat(paste(rep("\b", 70), collapse=""))

# Moment-matches a normal, truncated at x >= 0,
# with mean and variance m and v, respectively.
# following Phoebus Dhrymes'
# "Moments of Truncated (Normal) Distributions".
positive.moment.match = function(m, v) {
  m1 = -as.vector(m)
 v[v<0] = 1e10   # XXX hack
  s = sqrt(v)
  z = -m1 / s
  a = dnorm(z) / pnorm(z)

  # hack to deal with when z is very negative
  r = cbind(m = ifelse(z < -30, 1e-3, - (m1 - s * a)),
    v = ifelse(z < -30, 1e-6, v * (1 - z*a - a^2)))
#  r = cbind(m = - (m1 - s * a), v = v * (1 - z*a - a^2))
  r
}

# Messages constraining a variable to be positive
# (but using canonical parameterization.)
positive.moment.match.canonical = function(x) {
  mv = canonical.to.mean.and.variance(x)
  mm = positive.moment.match(mv[,"m"], mv[,"v"])
  mean.and.variance.to.canonical(mm)
}

# Original, slower, version of this.
mvnorm.2.diag = function(m, v, A, b, b.var) {
  V = Diagonal( x = v )

  # first, compute (A V A^T + v I) ^ -1, which I'll call B
  M = as.matrix( A %*% V %*% t(A) + Diagonal(x = b.var) )

  B = pseudoinverse(M)
#  B = chol2inv(chol(M))

  # ??? can we avoid computing the whole covariance here?
  cbind(m = m - as.vector(V %*% t(A) %*% B %*% (A %*% m - b)),
    v = as.vector(diag(V - V %*% t(A) %*% B %*% A %*% V)))
}

# Adds in a constraint of the form Ax ~ b, and returns the
# full covariance matrix.
mvnorm.2 = function(m, v, A, b, b.var) {
  V = Diagonal( x = v )

  # first, compute (A V A^T + v I) ^ -1, which I'll call B
  M = as.matrix( A %*% V %*% t(A) + Diagonal(x = b.var) )

  B = pseudoinverse(M)
#  B = chol2inv(chol(M))

  list(m = m - as.vector(V %*% t(A) %*% B %*% (A %*% m - b)),
    V = as.matrix( V - V %*% t(A) %*% B %*% A %*% V) )
}


# Distribution of x ~ N(m,v) | Ax ~ N(b,b.var).
lin.constraint.1 = function(m, v, A, b, b.var) {
# print("m =")
# print(m)
# print("v =")
# print(v)
  # first, compute (A V A^T + diag(v)) ^ -1, which I'll call B
#  M = as.matrix( A %*% (v * t(A)) + Diagonal(x = b.var) )
  M = A %*% (v * t(A))
  diag(M) = diag(M) + b.var
#  M = as.matrix( A %*% (v * t(A)) + Diagonal(x = b.var) )

  # using Cholesky would be slightly faster, but maybe less stable
  M.qr = qr(M)

# print("M =")
# print(M)
  r = cbind(m = m - v * as.vector(t(A) %*% solve(M.qr, A %*% m - b)),
    v = v - v * apply(A * solve(M.qr, A), 2, sum) * v)
# print(colnames(r))
  r
}

# Message constraining "Ax ~ N(b, b.var)".
lin.constraint.factor = function(A, b, b.var) function(x) {
  mv = canonical.to.mean.and.variance(x)
#  r = mvnorm.2.diag(mv[,"m"], mv[,"v"], A, b, b.var)
  r = lin.constraint.1(mv[,"m"], mv[,"v"], A, b, b.var)
# print("r =")
# print(r)
# print(colnames(r))
  mean.and.variance.to.canonical(r)
}

# Approximates a region constrained such that x >= 0.
# Args:
#   A, b, b.var - these give the constraint that
#     Ax ~ Normal(b, b.var), x >= 0
#     (For exact constraints, b.var can be 0.)
#   prior.var - variance for the prior (can be Inf)
# Returns: list with components
#   m, v - the mean and variance of the posterior
#   t - the prior times the terms (without the linear constraint)
#   update.stats - matrix with mean and variance of update sizes
approx.region = function(A, b, b.var, prior.var=Inf,
  converge.tolerance = 1e-9, max.iters=100) {
  n = ncol(A)

#  debug.dir = "~/tmp/approx.region.debug"
  debug.dir = NULL

  # prior (for now, restricted to be diagonal)
  prior = mean.and.variance.to.canonical(cbind(m=rep(0,n), v=rep(prior.var,n)))

  # the term approximations (initially flat)
  terms = mean.and.variance.to.canonical(cbind(m=rep(0,n), v=rep(Inf,n)))

  # the (soft) linear constraint
  lin.constraint = lin.constraint.factor(A, b, b.var)

  # the posterior
  q = mean.and.variance.to.canonical(cbind(m=rep(0,n), v=rep(1,n)))
#  q = lin.constraint( prior + terms )
#  q = prior

  # convergence statistics
  update.stats = NULL

  for(iter in 1:max.iters) {
    terms.old = terms

    terms.1 = q - terms

    if (!is.null(debug.dir)) {
      system(paste("mkdir -p", debug.dir))
      ep.trace = list(terms=terms, terms.1=terms.1, q=q)
      save(ep.trace, file=paste(debug.dir, "/", iter, ".Rdata", sep=""))
    }

    mm = positive.moment.match.canonical(terms.1)
#    mm[ is.nan(mm[,"e1"]) | is.nan(mm[,"e2"]) |
#      is.infinite(mm[,"e1"]) | is.infinite(mm[,"e2"]) , ] = 0

# print(as.vector(mm))
    # update terms.
#    terms = 0.5 * (mm - q) + 0.5 * terms   #  XXX not sure this is right
    terms = (mm - q) + terms

#    terms[ is.na(terms) ] = 0

    # add in Ax ~ N(-,-) constraint
    q = lin.constraint( prior + terms )

    # ??? show change in mean and variance separately?
    diff = apply(abs(canonical.to.mean.and.variance(terms.old) - canonical.to.mean.and.variance(terms)), 2, max)
cat(backspace, signif(diff, 2), " ")
    update.stats = rbind(update.stats, diff)

    # possibly stop early
    if (max(diff, na.rm=TRUE) <= converge.tolerance)
      break
  }

  mv = canonical.to.mean.and.variance(q)
  list(m = mv[,"m"], v = mv[,"v"],
    t = canonical.to.mean.and.variance( prior + terms ),
    update.stats = update.stats)
}

