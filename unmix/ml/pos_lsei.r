# Attempt at a simple greedy way of solving a particular
# lsei-type problem (in which the constraints are all "x >= 0".)

library(Matrix)
library(corpcor)
library(limSolve)

# Conditions a normal distribution (with diagonal covariance)
# on a linear constraint, somewhat as in the EP stuff.
# Args:
#   v - the (diagonal) variance
#   A, b.var - these define the "Ax = N(b, b.var)" constraint
#   b - also part of the constraint (separated out to avoid
#     recomputing the pseudoinverse
#   m - the current mean
# Returns: the updated mean
mvnorm.cond.diag = function(v, A, b.var) {
  V = Diagonal( x = v )

  # first, compute (A V A^T + v I) ^ -1, which I'll call B
  M = as.matrix( A %*% V %*% t(A) + Diagonal(x = b.var) )

  B = pseudoinverse(M)
#  B = chol2inv(chol(M))

  # we return a function which computes this
  function(m, b)
    m - as.vector(V %*% t(A) %*% B %*% (A %*% m - b))
}

# Goes "as far as possible" in some direction, without
# any element of a vector going negative.
move.bounded = function(a, b) {

  # the amount we would go in each direction
  delta = b - a

  s = - a / delta
  s[ is.na(s) ] = 0

  # If any element of s is negative, it doesn't matter, in
  # terms of hitting boundaries. However, if any element is
  # positive, we can only add that "amount" of delta.

  m = min(1, s[ s > 0 ])

  a + m * delta
}

# Finds a min-norm positive solution to Ax = b.
# Args:
#   A, b.var, b - these give the linear constriaint
#   max.iters - maximum number of iterations to go
#   converge.tol - if an update has a 2-norm smaller than
#     this, then we consider that we've "converged"
#     (even though the constraints may not be satisfied)
# Returns: smallest |x|, with all x >= 0, matching the
#   constraints (or the nearest non-negative point, if
#   there is no such point)
pos.lsei = function(A, b.var=NULL) {
  if (is.null(b.var)) {
    b.var = rep(0, nrow(A))
  }

  # function which takes a point to the nearest point
  # matching the constraints
  f = mvnorm.cond.diag(rep(1, ncol(A)), A, b.var)

  function(b, max.iters = 100, converge.tol = 1e-5) {


    # start with the "truncated pseudoinverse"
    x = pseudoinverse(A) %*% b
    x[ x < 0 ] = 0

    update.stats = cbind(err = max(abs(A %*% x - b)),
      amount.moved = NA)

    # matrix to hold the results
  #  X = matrix(NA, nrow=max.iters, ncol=length(x))
  #  X[1,] = x

    err = max(abs(A %*% x - b))

    for(i in 2:max.iters) {
cat(i, "")
      x1 = move.bounded(x, f(x, b))
      x1[x1 < 0] = 0
      err1 = max(abs(A %*% x1 - b))
      update.stats = rbind(update.stats,
        cbind(err = err1, amount.moved = max(abs(x1 - x))))
      x = x1

      if (abs(err1) <= converge.tol) {
 cat("converged in", i, "\n")
        return(x1)
      }
      x = x1
      err = err1
    }

    x
  }
}


