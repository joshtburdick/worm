# Approximates a constrained region using EP.

library(Matrix)
library(corpcor)

source("git/unmix/ept/normal.r")
# source("git/unmix/ept/matrix_inv_lemma.r")   XXX um, don't think I need this

backspace = paste(rep("\b", 70), collapse="")

# Moment-matches a normal, truncated at x >= 0,
# with mean and variance m and v, respectively.
# following Phoebus Dhrymes'
# "Moments of Truncated (Normal) Distributions".
positive.moment.match = function(m, v) {
  m1 = -as.vector(m)
#  v[v<0] = 1e10    # XXX hack
  s = sqrt(v)
  z = -m1 / s
  a = dnorm(z) / pnorm(z)

  # hack to deal with when z is very negative
  r = cbind(m = ifelse(z < -30, 1e-4, - (m1 - s * a)),
    v = ifelse(z < -30, 1e-10, v * (1 - z*a - a^2)))
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

#  cat("m in", range(m), "   v in", range(v), "\n")

  # first, compute (A V A^T + diag(v)) ^ -1, which I'll call B
#  M = as.matrix( A %*% (v * t(A)) + Diagonal(x = b.var) )
  M = A %*% (v * t(A))
  diag(M) = diag(M) + b.var
#  M = as.matrix( A %*% (v * t(A)) + Diagonal(x = b.var) )
# print(v[1:3])
# print("M =")
# print(M[1:3,1:3])
# print(as.numeric(determinant(M, logarithm=TRUE)))

  # using Cholesky would be slightly faster, but maybe less stable
r = NULL
if (FALSE) {  M.qr = qr(M)

  r = cbind(m = m - v * as.vector(t(A) %*% solve(M.qr, A %*% m - b)),
    v = v - v * apply(A * solve(M.qr, A), 2, sum) * v)
}

  M.pi = pseudoinverse(M)

  # trying to use the pseudoinverse here
  r = cbind(m = m - v * as.vector(t(A) %*% (M.pi %*% (A %*% m - b))),
    v = v - v * apply(A * (M.pi %*% A), 2, sum) * v)

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
    q.old = q
    terms.old = terms

    terms.1 = q - terms

    mm = positive.moment.match.canonical(terms.1)

# ???
#    mm = positive.moment.match.canonical(q) - terms

#    mm[ is.nan(mm[,"e1"]) | is.nan(mm[,"e2"]) |
#      is.infinite(mm[,"e1"]) | is.infinite(mm[,"e2"]) , ] = 0

# print(as.vector(mm))
    # update terms.
    # one way to add damping. XXX not sure this is right.
#    terms = 0.03 * (mm - q) + terms
    terms = (mm - q) + terms
    if (!is.null(debug.dir)) {
      system(paste("mkdir -p", debug.dir))
      ep.trace = list(terms=terms, terms.1=terms.1, q=q, mm=mm)
      save(ep.trace, file=paste(debug.dir, "/", iter, ".Rdata", sep=""))
    }

#    terms[ is.na(terms) ] = 0

    # add in Ax ~ N(-,-) constraint
    q = lin.constraint( prior + terms )

    # constraints on posterior
if (FALSE) {
# print(sum(is.na(q)))
    q1 = canonical.to.mean.and.variance(q)
    q1[ q1[,"m"] < 0 , "m" ] = 0
    q1[ q1[,"v"] < 1e-10, "v" ] = 1e-10
    q = mean.and.variance.to.canonical(q1)
# print(sum(is.na(q)))
}
    # ??? show change in mean and variance separately?
    diff = apply(abs(canonical.to.mean.and.variance(q) - canonical.to.mean.and.variance(q.old)), 2, max)
 cat(backspace, signif(diff, 2), " ")
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

# Also approximates a region constrained such that x >= 0,
# but includes damping (which may help convergence.)
# Args:
#   A, b - these give the constraint that
#     Ax = b, x >= 0
#   converge.tolerance - determines when to stop
#   max.iters - maximum number of iterations
#   prior.var - the prior variance (presumed diagonal)
#   damping - initial damping factor (no damping = 1)
#   damping.adjust - amount to multiply damping by,
#     if an update leads to NaN / Inf
# Returns: list with components
#   m, v - the mean and variance of the posterior
#   t - the prior times the terms (without the linear constraint)
#   update.stats - matrix with mean and variance of update sizes
approx.region.damping = function(A, b, converge.tolerance = 1e-9,
    max.iters=100, prior.var=Inf, damping=1, damping.adjust=0.5) {
  n = ncol(A)

  # prior (for now, restricted to be diagonal)
  prior = mean.and.variance.to.canonical(cbind(m=rep(0,n), v=rep(prior.var,n)))

  # the term approximations (initially flat)
  terms = mean.and.variance.to.canonical(cbind(m=rep(0,n), v=rep(Inf,n)))

  # the linear constraint (for now, this is exact)
  lin.constraint = lin.constraint.factor(A, b, 0 * b)

  # the posterior (??? initialize with pseudoinverse?)
  q = mean.and.variance.to.canonical(cbind(m=rep(0,n), v=rep(1,n)))
#  q = lin.constraint( prior + terms )

  # convergence statistics
  update.stats = NULL

  # definition of "numbers looking reasonably defined"
  num.defined = function(x)
    ((sum(is.na(x)) == 0) && (sum(is.nan(x)) == 0) && (sum(is.infinite(x)) == 0))

  for(iter in 1:max.iters) {
    error.flag = FALSE
    update.diff = NULL

    tryCatch({
      q.old = q
      terms.old = terms

      terms.1 = q - terms
      mm = positive.moment.match.canonical(terms.1)

      # update terms, with damping. XXX not sure this is right.
      terms.new = damping * (mm - q) + terms

      # add in Ax ~ N(-,-) constraint
      q.new = lin.constraint( prior + terms.new )

      # how much the posterior changed
      update.diff <- apply(abs(canonical.to.mean.and.variance(q.new) -
        canonical.to.mean.and.variance(q.old)), 2, max)

    }, error = function(e) error.flag <<- TRUE)

#    print(update.diff)
#    print(class(update.diff))
#    print(length(update.diff))
#    cat("m:", update.diff["m"], update.stats[nrow(update.stats),"m"], "\n")
#    cat("v:", update.diff["v"], update.stats[nrow(update.stats),"v"], "\n")

    # if that seemed to work, make the update
    if ((!error.flag) && num.defined(terms.new) && num.defined(q.new) &&
#        (min(canonical.to.mean.and.variance(q.new)[1,"m"]) >= 1e-6) &&
        (!is.null(update.stats)) &&
        (update.diff["m"] <= 2 * update.stats[nrow(update.stats),"m"]) &&
        (update.diff["v"] <= 2 * update.stats[nrow(update.stats),"v"])) {
      terms = terms.new
      q = q.new
    }
    else {
      if (iter >= 2)
        damping = damping * damping.adjust
    }

 cat(backspace, damping, update.diff, backspace)   # signif(diff, 2)
    update.stats = rbind(update.stats, c(damping=damping, update.diff))

    # possibly stop early
    if ((!error.flag) && (max(update.diff) <= converge.tolerance))
      break
  }
cat(" \n")

# print(update.stats)

  mv = canonical.to.mean.and.variance(q)
  list(m = mv[,"m"], v = mv[,"v"],
    t = canonical.to.mean.and.variance( terms ),
    update.stats = update.stats)
}


