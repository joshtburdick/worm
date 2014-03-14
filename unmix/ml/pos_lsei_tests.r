# Tests of pos_lsei functions.

source("git/unmix/ml/pos_lsei.r")



r0 = pos.lsei(t(c(1,1,1)))(1)

set.seed(12)
# A1 = matrix(runif(50000), nrow=50)
A1 = matrix(rbeta(50000, 0.1, 0.1), nrow=50)
x1 = rgamma(1000, shape=0.1, scale=50)
# x1 = rbeta(1000, 0.1, 0.1)
b1 = A1 %*% x1

print(system.time(r1 <- pos.lsei(A1)(b1)))

err = (A1 %*% r1) - as.vector(b1)

sse = apply(err^2, 2, max)

# error in meeting the constraints
print(sse)

print(system.time(r2 <- ldei(A1, b1, diag(ncol(A1)), rep(0, ncol(A1)))))
# print(system.time(r2 <- nnls(A1, b1)))

# difference from the ldei() solution
print(sum( (r1 - r2$X)^2 )) 


