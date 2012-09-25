
library(corpcor)

source("git/unmix/ept/dirichlet.r")








# Constrains one variable to be a linear function of another.
lin.factor = function(A, b) {
  Ap = pseudoinverse(A)
  function(m) {
    x1 = canonical.to.mean.and.variance(m$x)
    y1 = canonical.to.mean.and.variance(m$y)
    m.y = cbind(m = A %*% x1[,"m"] + b, v = t(A) %*% diag(x1[,"v"]) %*% A)
    m.x = cbind(m = Ap %*% (y1[,"m"] - b), v = t(Ap) %*% diag(y1[,"v"]) %*% Ap)
    list(m = list(x = mean.and.variance.to.canonical(m.x),
      y = mean.and.variance.to.canonical(m.y)), le=0)
  }
}




