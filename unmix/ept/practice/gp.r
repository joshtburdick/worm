# Gaussian Processes.

library(Matrix)

# Kernels all take as input two matrices, each of which has
# one example per row.

# Gaussian RBF kernel.
gaussian.rbf.kernel = function(length) function(x1, x2) {
  x = rbind(x1, x2)
  d1 = dim(x1)[1]
  d2 = dim(x2)[1]
  k = exp( -0.5 * ( as.matrix(dist(x))^2 / length ) )
  k[1:d1,(d1+1):(d1+d2)]
}

# Polynomial kernel.
# ??? does this need converting to a distance to make
# it SPD?
poly.kernel = function(a, b, k) function(x1, x2)
  (a + b * (x1 %*% t(x2))) ^ k

# Log-likelihood, for some data set.
# ??? does the "full Bayesian" estimate of this relate
# to cross-validation?
gp.ll.dense = function(x, y, kernel.fun, si2) {
  k = kernel.fun(x, x) + diag(si2, dim(x)[1])
  u = chol(k)
  z = backsolve(u, backsolve(t(u), y, upper.tri=FALSE))
  n = length(y)

  - sum(log(diag(u))) - 0.5 * sum( y * z ) - (n/2)*log(2*pi)
}

# Log-likelihood, for some data set, allowing sparsity.
# We assume that the diagonal of 
gp.ll = function(kernel.matrix, y, si2) {
  n = length(y)
  k = kernel.matrix + diag(si2, n)
  chol(k)

  -0.5 * (determinant(k, logarithm=TRUE)$modulus[1] +
    sum(y * solve(k, y)) +
    n * log(2*pi))
}

# Gaussian process posterior predictive distribution,
# approximated a la MacKay and _____.
# gp.pp.approx = function(



# Gaussian process posterior predictive distribution.
# Doesn't use sparse matrices; not intended to be fast.
# Args:
#   x - predictors
#   y - what's being predicted (currently a vector;
#     should have as many entries as x has rows)
#   kernel.fun - function which returns a kernel matrix
#   si2 - "noise" variance
# Returns:
#   matrix with two columns
gp.pp.1 = function(x, y, kernel.fun, si2) {
  k = kernel.fun(x, x) + diag(si2, dim(x)[1])
  u = chol(k)
#  y1 = backsolve(u, backsolve(t(u), y, upper.tri=FALSE))
  y1 = solve(k, y)

  function(x1) {
print(dim(k))
    k1 = as.matrix(kernel.fun(x1, x))
print(dim(k1))
#    z = backsolve(u, backsolve(t(u), k1, upper.tri=FALSE))

    # compute full covariance
    si2.full = kernel.fun(x1,x1) - k1 %*% solve(k, t(k1))

    list(mu = k1 %*% y1, si2 = diag(si2.full))
#    list(mu = t(k1) %*% y1)  #, u=u, k=k, k1=k1, y1=y1, z=z)
#      si2 = kernel.fun(x1,x1) - t(k1) %*% z )
#   ??? as.matrix(x1)
  }
}

# Some sample matrices for testing.
set.seed(41)
x1a = matrix(runif(9), nrow=3)
x1b = matrix(runif(6), nrow=2)

set.seed(42)
x2 = matrix(rnorm(1000), nrow=100)
y2 = x2 %*% runif(10) + rnorm(100, sd=0.1)

set.seed(42)
x3 = as.matrix( c(1:100)/10 )
y3 = as.matrix( sin(x3) + 0.5*(x3/8)^2.5 - log(x3/2) + rnorm(100, sd=0.3) )

# add an outlier by hand
y3[10] = y3[10] - 3

# Plots some data points, the regressed points, and 
# the mean and 2 s.d. bounds of the estimate.
# Args:
#   length.scale - scale for the Gaussian RBF kernel
#   noise.si2 - the noise variance
test1 = function(length.scale=1, noise.si2=1) {
  r = gp.pp.1(x3,y3,gaussian.rbf.kernel(length.scale), noise.si2)(x3)
#  foo = gp.pp(x3, y3, poly.kernel(1, 1, 5), noise.si2)(x3)
  xlim = range(x3)
  ylim = range(y3, r$mu)
  plot(x3, y3, xlim=xlim, ylim=ylim, type="p")
  par(new=TRUE)
  plot(x3, r$mu - 2 * sqrt(r$si2),
    xlim=xlim, ylim=ylim, xlab="", ylab="",
    type="l", col="#b0b0b0")
  par(new=TRUE)
  plot(x3, r$mu + 2 * sqrt(r$si2),
    xlim=xlim, ylim=ylim, xlab="", ylab="",
    type="l", col="#b0b0b0")
  par(new=TRUE)
  plot(x3, r$mu, xlim=xlim, ylim=ylim, xlab="", ylab="",
    type="l", col="#808080")
  r
}

# a simple test
if (FALSE) {
  pdf("git/unmix/ept/practice/gp.pdf")
  foo = test1(0.2, 0.2)
  dev.off()
}

