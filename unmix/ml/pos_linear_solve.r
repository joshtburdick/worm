# Solves many constrained linear systems at once.

library(corpcor)

source("git/utils.r")

# Finds the closest solution to AX = B.
# Args:
#   A, B - these give the constraint
# Returns: function which finds the nearest X (note
# that this avoids recomputing A's pseudoinverse)
lin.constraint = function(A, B) {
  P = t(A) %*% pseudoinverse( A %*% t(A) )
  function(X) {
    X - P %*% ( A %*% X - t(B) ) 
  }
}

# Goes "as far as possible" from A to B, without
# any element of a vector going negative.
# (Each column is treated separately.)
move.pos = function(A, B) {
  D = B - A
  S = - A / D
  S[ is.na(S) ] = 0

  # for each column of S, find the smallest positive
  # entry (that's how far we could "go")
  S[ S < 0 ] = Inf
  m = apply(S, 2, min)

  # we don't want to "overshoot" B
  m = pmin(1, m)

  # "go" that amount, in the direction towards B
  # XXX the transposes here may be slow
  A + t( m * t(D) )
}

# Solves for X in AX = B, X >= 0 (in other words, solves
# for many variables at once.)
# Args:
#   A, B - these give the matrix
#   X - the starting point
#   max.iters, eps - stopping criteria
# Returns: a positive solution of X (or, if there isn't such,
#   then the closest positive solution by least squares.)
pos.linear.solve = function(A, B, X, max.iters=50, eps=1e-10) {
  update.stats = NULL
  lc = lin.constraint(A, B)

  for(iter in 1:max.iters) {

    # move toward the closest point satisfying the constraints
    X1 = move.pos(X, lc(X))

    update.size = max(abs(X1 - X), na.rm=TRUE)
write.status(paste(iter, update.size, "\n"))
    update.stats = rbind(update.stats,
      data.frame(update.size = update.size))
    X = X1

    if (update.size <= eps)
      return( list(X=X, update.stats=update.stats) )
  }

  list(X = X, update.stats = update.stats)
}




