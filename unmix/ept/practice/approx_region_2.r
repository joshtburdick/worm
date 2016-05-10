# Approximates a constrained region using EP.

library(Matrix)
library(corpcor)

source("git/unmix/ept/normal.r")

source("git/unmix/ept/gamma.r")

# source("git/unmix/ept/matrix_inv_lemma.r")   XXX um, don't think I need this

# trying this variant of the conditional thing ?
# doesn't seem to work
source("git/unmix/ept/practice/gamma_conditional_7.r")

backspace = paste(rep("\b", 70), collapse="")

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
  r = cbind(m = ifelse(z < -30, 0, - (m1 - s * a)),
    v = ifelse(z < -30, 0, v * (1 - z*a - a^2)))
#  r = cbind(m = - (m1 - s * a), v = v * (1 - z*a - a^2))
  r
}

# Messages constraining a variable to be positive
# (but using canonical gamma parameterization.)
positive.moment.match.canonical.gamma = function(x) {
  mv = gamma.n2mv(t(x))
  mm = positive.moment.match(mv["m",], mv["v",])
#print(mm)
#  i = is.na(mm[,"v"])
  r = t(gamma.mv2n(t(mm)))
#  r[i,"e1"] = -1
#  r[i,"e2"] = 0
#print(r)
  r
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

# Distribution of one of these, conditional on the others (hack).
mvnorm.cond.other = function(m, v, A, b) {
  n = length(m)
  r = list(m = rep(NA, n), v = rep(NA, n))

  for(i in 1:n) {
    m1 = m
    v1 = v
    m1[i] = 0
    v1[i] = 1e10
    a = mvnorm.2.diag(m1, v1, A, b, 0*b)
    r$m[i] = a[i,"m"]
    r$v[i] = a[i,"v"]
  }
  r
}

# Like the above, but also moment-matching to be > 0.
mvnorm.cond.other.pos = function(x, A, b) {
  r = mvnorm.cond.other(x["m",], x["v",], A, b)

  t(positive.moment.match(r$m, r$v))
}

mvnorm.cond.other.pos.gamma = function(A, b)
  function(x)
    gamma.mv2n(mvnorm.cond.other.pos(gamma.n2mv(x), A, b))





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

# Linear constraint, with conversions to/from gamma params
# (but no moment matching.)
lin.constraint.gamma = function(A, b, b.var) function(x) {
  mv = gamma.n2mv(x)
# print(mv)
#  r = mvnorm.2.diag(mv[,"m"], mv[,"v"], A, b, b.var)
  r = lin.constraint.1(mv["m",], mv["v",], A, b, b.var)
#  print("r =")
#  print(r)
  gamma.mv2n(t(r))
}

# Linear constraint, with conversions to/from gamma params,
# and moment matching to the positive region.
lin.constraint.gamma.pos = function(A, b, b.var) function(x) {
  mv = gamma.n2mv(x)
  r = t(lin.constraint.1(mv["m",], mv["v",], A, b, b.var))
  print("r =")
print(r)
  mm = t(positive.moment.match(r["m",], r["v",]))
print(mm)
  gamma.mv2n(mm)
}

# Approximates a region constrained such that x >= 0.
# Uses the Gamma distribution internally.
# Args:
#   A, b, b.var - these give the constraint that
#     Ax ~ Normal(b, b.var), x >= 0
#     (For exact constraints, b.var can be 0.)
#   prior.var - variance for the prior (can be Inf)
# Returns: list with components
#   m, v - the mean and variance of the posterior
#   t - the prior times the terms (without the linear constraint)
#   update.stats - matrix with mean and variance of update sizes
approx.region.gamma = function(A, b, b.var, prior.var=Inf,
  converge.tolerance = 1e-9, max.iters=100) {
  n = ncol(A)

#  debug.dir = "~/tmp/approx.region.debug"
  debug.dir = NULL

  # prior (for now, restricted to be diagonal)
  prior = 0 * t(gamma.mv2n(rbind(m=rep(1,n), v=rep(1,n))))

  # the term approximations (initially flat?)
#  terms = 0 * prior
  terms = t(gamma.mv2n(rbind(m=rep(1,n), v=rep(1,n))))

  # the (soft) linear constraint
  lin.constraint = lin.constraint.gamma(A, b, b.var)

  # the posterior
#  q = mean.and.variance.to.canonical(cbind(m=rep(0,n), v=rep(1,n)))
  q = lin.constraint( prior + terms )
#  q = prior

  # convergence statistics
  update.stats = NULL

  for(iter in 1:max.iters) {
# print(iter)
# print(t(gamma.n2mv(t(q))))

    terms.old = terms

    terms.1 = q - terms
# print("terms.1")
# print(terms.1)
    mm = positive.moment.match.canonical.gamma(terms.1)
  mm[is.na(mm)] = 0
print("mm")
print(mm)
#    mm = positive.moment.match.canonical.gamma(q) - terms
#    mm = terms.1

#    mm[ is.nan(mm[,"e1"]) | is.nan(mm[,"e2"]) |
#      is.infinite(mm[,"e1"]) | is.infinite(mm[,"e2"]) , ] = 0

# print(as.vector(mm))
    # update terms.
    # one way to add damping. XXX not sure this is right.
#    terms = 0.5 * (mm - q) + 0.5 * terms
    terms = (mm - q) + terms
 print("terms")
 print(terms)

    if (!is.null(debug.dir)) {
      system(paste("mkdir -p", debug.dir))
      ep.trace = list(terms=terms, terms.1=terms.1, q=q, mm=mm)
      save(ep.trace, file=paste(debug.dir, "/", iter, ".Rdata", sep=""))
    }

#    terms[ is.na(terms) ] = 0

    # add in Ax ~ N(-,-) constraint
    q = positive.moment.match.canonical.gamma( lin.constraint( prior + terms ))
#    q = lin.constraint(prior+terms)
#    q = prior+terms

    # constraints on posterior
if (FALSE) {
print(sum(is.na(q)))
    q1 = canonical.to.mean.and.variance(q)
    q1[ q1[,"m"] < 0 , "m" ] = 0
    q1[ q1[,"v"] < 1e-10, "v" ] = 1e-10
    q = mean.and.variance.to.canonical(q1)
print(sum(is.na(q)))
}
    # ??? show change in mean and variance separately?
    diff = apply(abs(canonical.to.mean.and.variance(terms.old) - canonical.to.mean.and.variance(terms)), 2, max)
cat(backspace, signif(diff, 2), " ")
    update.stats = rbind(update.stats, diff)

    # possibly stop early
    if (max(diff, na.rm=TRUE) <= converge.tolerance)
      break
  }

#  mv = t(gamma.n2mv(t(q)))
# ???
  mv = t(gamma.n2mv( t(lin.constraint.gamma(A, b, b.var)(q) )))

  list(m = mv[,"m"], v = mv[,"v"],
    t = gamma.n2mv( t(prior + terms) ),
    update.stats = update.stats)
}

# XXX hack
pmm = function(x) t(positive.moment.match.canonical.gamma(t(x)))

# Another take at the same thing.
# (The variance seems wrong for this, usually).
approx.region.gamma.2 = function(A, b, b.var,
    converge.tolerance=1e-9, max.iters=100) {
  n = ncol(A)

  # the "marginal" terms
  t.m = gamma.mv2n(rbind(m=rep(1,n), v=rep(1,n)))

  # the "constraint" term
#  t.c = gamma.mv2n(rbind(m=rep(1,n), v=rep(1,n)))

  # the linear constraints
  lc = lin.constraint.gamma(A, b, b.var)
  lcp = lin.constraint.gamma.pos(A, b, b.var)

  # the posterior
  q = lcp(t.m)

  # convergence statistics
  update.stats = NULL

  for(iter in 1:max.iters) {
#  print(gamma.n2mv(terms))

    t.m = q - lc(t.m)   # lc(t.m)

    # update posterior: terms, with the Ax ~ N(-,-) constraint
    q.new = lcp(t.m)
# print(q.new)

    # FIXME this should be based on change in the terms, properly
    diff = apply(abs(gamma.n2mv(q.new) - gamma.n2mv(q)), 1, max)
cat(backspace, signif(diff, 2), " ")
    update.stats = rbind(update.stats, diff) 

    # possibly stop early
    if (max(diff, na.rm=TRUE) <= converge.tolerance)
      break

    q = q.new
# print(gamma.n2mv(q))
  }
cat("\n")

  mv = gamma.n2mv(lc(q))

  list(m = mv["m",], v = mv["v",],
    t = gamma.n2mv( t.m ),
    update.stats = update.stats)
}

# Yet another take on the same thing.
approx.region.gamma.3 = function(A, b, b.var,
    converge.tolerance=1e-9, max.iters=100) {
  n = ncol(A)

  # the "marginal" terms
#  t.m = gamma.mv2n(rbind(m=rep(1,n), v=rep(1,n)))

  # the "constraint" term
  t.c = gamma.mv2n(rbind(m=rep(1,n), v=rep(1,n)))

  # the linear constraints
  lc = lin.constraint.gamma(A, b, b.var)
  lcp = lin.constraint.gamma.pos(A, b, b.var)

  lcg = gamma.cond.1(A, b)

  # the posterior
#  q = t.m
#  t.c = q - t.m
  q = t.c
  t.c = 0 * t.c

  # convergence statistics
  update.stats = NULL

  for(iter in 1:max.iters) {
#  print(gamma.n2mv(terms))
print(iter)
    t.m = q - t.c  # was t.m

    # update posterior: terms, with the Ax ~ N(-,-) constraint
#    q.new = lcp(q - t.c)
# XXX trying new constraint code
     q.new = lcg(q - t.c)

    t.c = t.c + (q.new - q)
# print(q.new)

    # FIXME this should be based on change in the terms, properly
    diff = apply(abs(gamma.n2mv(q.new) - gamma.n2mv(q)), 1, max)
cat(backspace, signif(diff, 2), " ")
    update.stats = rbind(update.stats, diff) 

    # possibly stop early
    if (max(diff, na.rm=TRUE) <= converge.tolerance)
      break

    q = q.new
# print(gamma.n2mv(q))
  }
cat("\n")

  mv = gamma.n2mv(q)

  list(m = mv["m",], v = mv["v",],
    t = gamma.n2mv( t.c ),
    update.stats = update.stats)
}

if (TRUE) {

# XXX for practice
n = 100
A = t(rep(1,n))
t.m = gamma.mv2n(rbind(m=rep(1,n), v=rep(1,n)))

lc = lin.constraint.gamma(A, 1, 0)

r = approx.region.gamma.3(A, c(1), c(0))
}

if (TRUE) {
A1 = matrix(c(1,0.1,0.1,1,1,1), nrow=2)
r1 = approx.region.gamma.3(A1, c(1,3), c(0,0))

}

