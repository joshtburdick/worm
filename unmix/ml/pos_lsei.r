# Attempt at a simple greedy way of solving a particular
# lsei-type problem (in which the constraints are all "x >= 0".)

library(corpcor)
library(limSolve)

# Efficiency be darned! Using slower, old, mvnorm.2.diag().
source("git/unmix/ept/approx_region.r")






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

# ??? does this "get stuck" after one step?
# Apparently not.
find.mean.greedy = function(A, b, num.iters=10) {
  P = pseudoinverse(A)

  # start with the "truncated pseudoinverse"
  x = P %*% b
  x[ x < 0 ] = 0

  # matrix to hold the results
  X = matrix(NA, nrow=num.iters, ncol=length(x))
  X[1,] = x
  for(i in 2:num.iters) {
    r = mvnorm.2.diag(x, 1+0*x, A, b, 0*b)
    x.p = r[,1]
    x = move.bounded(x, x.p)
    x[x < 0] = 0
    X[i,] = x
  }

  X
}

