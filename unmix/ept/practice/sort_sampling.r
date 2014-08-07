# Trying to sample the sort matrix, and expression.
# Efficiency not a priority here.
# FIXME: add constraint so that rows of x add up to
# (or average) 1?

library("MASS")

# useful for printing progress
write.status = function(s) {
  max.length = 70
  s = substr(s, 1, max.length)

  bs = paste(rep("\b", max.length), collapse="")
  pad = paste(rep(" ", max.length - nchar(s)), collapse="")

  cat(paste0(bs, s, pad))
}

# Metropolis update rule.
metropolis = function(x, x1, ll) {
  a = ll(x)
  b = ll(x1)
  (b >= a) || rbinom(1, 1, exp(b-a))
}

# Unnormalized log-Dirichlet density.
ld.dirichlet = function(alpha, x) {
  sum((alpha - 1) * log(x))
}

# Log-likelihood for a matrix with rows
# from a Dirichlet.
dirichlet.ll = function(alpha) function(A) {
  s = 0
  for(i in 1:nrow(alpha))
    s = s + ld.dirichlet(alpha[i,], A[i,])
  s
}

# Finds the closest point in a hyperplane to a
# given point (essentially as seen in the EP code)
# Args:
#   A, b - these give the constraint
#   x - the point in question
# Returns: closest point to x in "Ax=b" hyperplane
closest.point = function(A, b, x) {
  M = chol2inv(chol(A %*% t(A)))
  x - as.vector(t(A) %*% M %*% (A %*% x - b))
}

# Utility to jump a random amount in some direction.
random.jump = function(x, d) {
  r = x / d
  lo = max(c(-Inf, -r[r>0]))
  hi = min(c(Inf, -r[r<0]))

  x + runif(1, lo, hi) * d
}

# Coordinate-directions sampling rule.
cda.sample = function(Z, x, num.iters) {
  for(i in 1:num.iters) {
    z = Z[ , sample(ncol(Z), 1) ]
    x = random.jump(x, z)
  }

  x
}

# Updates A and x.
# XXX if this "gets stuck", we simply don't update them.
# This is incorrect, but may be useful for debugging.
update.A.and.x = function(A, x, B, sd) {
  iter = 0
  while (iter < 100) {
    iter = iter + 1

    # sample A
    A1 = A * exp(rnorm(nrow(A) * ncol(A), mean=0, sd=sd))
    A1 = A1 / apply(A1, 1, sum)

    # update x to the "closest point to x in that hyperplane"
    x1 = x
    for(j in 1:ncol(x)) {
      x1[,j] = closest.point(A1, B[,j], x[,j])
    }

    # if all x's are positive, return the updated things
    if (all(x1 >= 0)) {
      return(list(A = A1, x = x1, iter = iter))
    }
  }

  # XXX not correct
  return(list(A = A, x = x, iter = iter))
}

# Samples expression, and the sort matrix.
# Args:
#   alpha - Dirichlet prior for rows of A
#   B - the expression data (whose rows should
#     sum to 1)
#   num.samples - number of samples to return
# Returns: list with components
#   A - samples of A
#   x - samples of x
sort.sample.1 = function(alpha, B, A, x, num.iters) {

  A.samples = array(dim=c(nrow(alpha), ncol(alpha), num.iters))
  x.samples = array(dim=c(nrow(x), ncol(x), num.iters))
  A.update.iters = rep(NA, num.iters)

  for(iter in 1:num.iters) {
write.status(iter)
    # save this sample
    A.samples[,,iter] = A
    x.samples[,,iter] = x

    # proposed update to A (and x)
    r = update.A.and.x(A, x, B, 0.1)
    A = r$A
    x = r$x
    A.update.iters[iter] = r$iter

    # Metropolis step, incorporating the prior
    if (metropolis(A, r$A, dirichlet.ll(alpha))) {
      A = r$A
      x = r$x
    }

    # update x in the new plane
if (FALSE) {
    Z = Null(t(A))
    for(j in 1:ncol(x)) {
      x[,j] = cda.sample(Z, x[,j], 100)
    }
}
  }

  list(A = A.samples, x = x.samples, A.update.iters = A.update.iters)
}

# Version of the above which picks a starting point.
sort.sample = function(A.params, B, num.iters) {

  # pick starting points
  A0 = matrix(rgamma(nrow(A.params)*ncol(A.params), 1, 1),
    nrow=nrow(A.params), ncol=ncol(A.params))
  A0 = A0 / apply(A0, 1, sum)
  x0 = matrix(nrow=ncol(A), ncol=ncol(B))

  for(j in 1:ncol(B)) {
    x0[,j] = ldei(E=A0, F=B[,j], G=diag(ncol(A0)), H=rep(0, ncol(A0)))$X
  }

  sort.sample.1(A.params, B, A0, x0, num.iters)
}
