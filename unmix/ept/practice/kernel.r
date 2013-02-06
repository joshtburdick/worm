# Various kernels, including fast-ish approximations.

library("Matrix")
library("FNN")

# Gaussian RBF kernel.
# ??? call this "squared-exponential"?
gaussian.rbf.kernel = function(length) function(x1, x2) {
  x = rbind(x1, x2)
  d1 = dim(x1)[1]
  d2 = dim(x2)[1]
  k = exp( -0.5 * ( as.matrix(dist(x))^2 / length ) )
  k[1:d1,(d1+1):(d1+d2)]
}

# Constructs a nearest neighbor matrix, using just the
# nearest neighbors for each data point.
# This version of it includes computes the whole square
dist.to.gaussian.rbf.kernel = function(knn.dist, k, length.scale) {
  num.rows = dim(knn.dist$nn.index)[1]

  # XXX if we re-do this for "distances between two sets of things",
  # this will need changing
  num.columns = num.rows

  # this tweak is needed so that the matrix is guaranteed symmetric
  i = rep(1:num.rows, k)
  j = as.vector(knn.dist$nn.index[,1:k])
  x = exp( -0.5 * as.vector(knn.dist$nn.dist[,1:k])^2 / length.scale)
  a = unique(data.frame(i = c(i,j), j = c(j,i), x = c(x,x)))

  m = sparseMatrix(dims = c(num.rows, num.columns), i = a$i, j = a$j, x = a$x)
  diag(m) = 1
  m = as(m, "symmetricMatrix")
}




