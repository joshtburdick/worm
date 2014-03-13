# Tests of pos_lsei functions.

source("git/unmix/ml/pos_lsei.r")



r0 = find.mean.greedy(t(c(1,1,1)), 1)

set.seed(423)
A1 = matrix(runif(50000), nrow=50)
x1 = rgamma(1000, shape=0.5, 1)
b1 = A1 %*% x1

print(system.time(r1 <- find.mean.greedy(A1, b1, num.iters=10)))

print(system.time(r2 <- ldei(A1, b1, diag(ncol(A1)), rep(0, ncol(A1)))))

err = (A1 %*% t(r1)) - as.vector(b1)

sse = apply(err^2, 2, sum)

# error in meeting the constraints
print(sse)

# difference from the ldei() solution
print(sum( (r1[nrow(r1),] - r2$X)^2 )) 


