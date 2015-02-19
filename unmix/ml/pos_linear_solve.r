# Solves many constrained linear systems at once.

library(corpcor)

source("git/utils.r")

# Finds the closest solution to AX = B.
# Args:
#   A, B - these give the constraint
# Returns: function which finds the nearest X (note
# that this avoids recomputing A's pseudoinverse)
lin.constraint = function(A, B) {
# was:
#  P = t(A) %*% pseudoinverse( A %*% t(A) )
  Ap = pseudoinverse(A)
  P = t(A) %*% t(Ap) %*% Ap
  function(X) {
    X - P %*% ( A %*% X - B ) 
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
#   X - the starting point (if NULL, the "truncated pseudoinverse"
#     is used)
#   normalize - indicates whether to normalize "rows", "columns"
#     or (the default, if "none") neither
#   max.iters, eps - stopping criteria
# Returns: a positive solution of X (or, if there isn't such,
#   then the closest positive solution by least squares.)
pos.linear.solve = function(A, B, X=NULL, normalize="none",
    max.iters=50, eps=1e-13) {
  update.stats = NULL
  lc = lin.constraint(A, B)

  # if no "starting point" is given, use a point "near zero"
  # ??? or "truncated pseudoinverse"?
  if (is.null(X)) {
    X = lc(matrix(0, nrow=ncol(A), ncol=ncol(B)))
    X[ X < 0 ] = 0
  }

  for(iter in 1:max.iters) {

    # move toward the closest point satisfying the constraints
    X1 = move.pos(X, lc(X))

    # also possibly normalize rows or columns
    if (normalize == "rows") {
      X1 = X1 / apply(X1, 1, sum)
      X1[ X1 < 0 ] = 0
    }
    if (normalize == "columns") {
      X1 = t( t(X1) / apply(X1, 2, sum) )
      X1[ X1 < 0 ] = 0
    }

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

