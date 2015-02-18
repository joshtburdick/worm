# Some tests of pos_linear_solve.r .

source("git/unmix/ml/pos_linear_solve.r")

library("limSolve")

# Very small smoke test of this.
# Args:
#   A, B - these give the problem
# Returns: the error, and minimum value
test1 = function(A, B) {

  X = pos.linear.solve(A, B)$X

  c(err = max(abs(A %*% X - B)), min=min(X))
}

# Comparison with lsei().
lsei.test.1 = function(A, B) {
  A = A / apply(A, 1, sum)
  B = B / apply(B, 1, sum)

  X = pos.linear.solve(A, B)$X

  X.lsei = X + NA    # XXX
  for(j in 1:ncol(B)) {
    r = ldei(A, B[,j], diag(ncol(A)), rep(0, ncol(A)))
    X.lsei[,j] = r$X
  }

  c(err = max(abs(A %*% X - B)), min=min(X),
    ldei.diff = max(abs(X - X.lsei)))
}

A1 = matrix(rgamma(40, shape=1, rate=1), nrow=4)
X1 = matrix(rgamma(30, shape=1, rate=1), nrow=10)
B1 = A1 %*% X1
lsei.test.1(A1, B1)

A2 = matrix(rgamma(400, shape=1, rate=1), nrow=4)
X2 = matrix(rgamma(300, shape=1, rate=100), nrow=100)
B2 = A2 %*% X2
print(lsei.test.1(A2, B2))

for(iter in 1:3) {
  A3 = matrix(rgamma(20000, shape=1, rate=1), nrow=20)
  X3 = matrix(rgamma(30000, shape=1, rate=100), nrow=1000)
  B3 = A3 %*% X3
  print(lsei.test.1(A3, B3))
}

# can it do least squares? what happens if the system
# is _over_determined?
A4 = matrix(rgamma(10, shape=1, scale=1), nrow=5)
X4 = matrix(rgamma(10, shape=1, scale=1), nrow=2)
B4 = A4 %*% X4
r4 = pos.linear.solve(A4, B4)


