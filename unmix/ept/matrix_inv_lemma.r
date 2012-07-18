# Inverts certain matrices fast, using the
# Sherman-Morrison-Woodbury lemma.

library(Matrix)

# Computes the diagonal of the inverse of A + X B Xt
# Args:
#   A.diag - diagonal of A
#   X - a matrix with many more rows than columns
#   B - a s.p.d. matrix
# Returns: diagonal of inverse of A + X B Xt
matrix.inv.diag.narrow = function(A.diag, X, B) {

  # inverse of A
  A.inv = Diagonal( x = 1 / A.diag )
  A.inv.diag = 1 / A.diag

  # intermediate result 
  C = chol2inv(chol(B)) + t(X) %*% A.inv %*% X

  # ??? can we avoid computing this whole matrix?
#  diag( A.inv - A.inv %*% X %*% solve(C, t(X)) %*% A.inv )
  A.inv.diag - A.inv.diag * apply(X * t(solve(C, t(X))), 1, sum) * A.inv.diag
}

if (FALSE) {
# Tests the above, by comparing with the naive function.
# Args: A, X, and B, as for matrix.inv.diag.narrow()
# Returns: mean-squared difference between the above
#   function's output, and the same thing computed naively
matrix.inv.diag.narrow.test = function(A.diag, X, B) {
  r = matrix.inv.diag.narrow(A.diag, X, B)

  r.naive = diag(chol2inv(chol( diag(A.diag) + X %*% B %*% t(X) )))

  sum( (r - r.naive)^2 )
}

set.seed(0)
t1 = list(A.diag = sample(10, 5),
  X = matrix(runif(10), nrow=5, ncol=2),
  B = matrix(c(2, -0.3, -0.3, 2), nrow=2))




}

