# Another take on this; somewhat of a smoke test.

source("git/unmix/ept/gamma.r")


# Rescales things to sum to something.
gamma.rescale.sum = function(x, b) {
  x1 = gamma.n2mv(x)
  
  m1 = sum(x1["m",])

  # rescale
  s = b / m1
  x1["m",] = s * x1["m",]
  x1["v",] = s * s * x1["v",]

  gamma.mv2n(x1)
}

# Flat distribution (normalized to have expected sum of 1,
# which shouldn't really matter).
gamma.flat = function(n) {
  x = gamma.mv2n(rbind(m=rep(1,n), v=rep(1,n)))
  gamma.rescale.sum(x, 1)
}

# Scales proportions of gamma-distributed variables.
# Args:
#   x - gamma distributions to scale
#   A - boolean giving a subset of x (XXX currently
#     this must include, and exclude, at least 2 vars)
#   b - proportion to give the "A == TRUE" set
# Returns: those gamma parameters, rescaled
gamma.cond.scale.1 = function(x, A, b) {
  # ??? add "drop=FALSE"?
  x[,A] = gamma.rescale.sum(x[,A], b)
  x[,!A] = gamma.rescale.sum(x[,!A], 1-b)
  x
}

# Yet another approximation of Ax = b, x > 0.
approx.region.gamma.1 = function(A, b) {

  # the approximating distributions
  p = array(0, dim=c(2, ncol(A), nrow(A)))
  dimnames(p)[[1]] = c("e1", "e2")

  # the prior
  q0.shape = 1
  q0.s = rbind(a = rep(q0.shape, ncol(A)), b = rep(1, ncol(A)))
  q0 = gamma.rescale.sum(gamma.s2n(q0.s), 1)

  # the posterior
  q = q0

  for(iter in 1:50) {

    # loop through the factors
    for(i in 1:nrow(A)) {

      # subtract out this term
      p1 = q - p[,,i]

      # add in this constraint
      p2 = gamma.cond.scale.1(p1, A[i,], b[i])
      
      # update this term
      p[,,i] = p1 + 0.01 * (p2 - q)

      # update posterior
      q = apply(p, c(1,2), sum)

    }


print(iter)
print(gamma.n2mv(q))  # XXX for debugging
  }

  q
}

A0 = rbind(c(T,T,T,F,F), c(F,F,T,T,T))
b0 = c(0.3, 0.8)
r0 = approx.region.gamma.1(A0, b0)

if (FALSE) {
set.seed(1)
A1 = rbind(1:10 >= 7, 1:10 %% 3 != 0)
b1 = c(0.5, 0.5)
r1 = approx.region.gamma.1(A1, b1)

cat("\n\n")
}
# this runs, and ends up with mean exactly equal to s.d.,
# which seems incorrect in general

if (FALSE) {
set.seed(42)
A2 = matrix(runif(100) >= 0.5, ncol=20)
x2 = runif(20)
x2 = x2 / sum(x2)
b2 = A2 %*% x2

r2 = approx.region.gamma.1(A2, b2)
# this doesn't converge
}

# trying to create an example which has sd != sqrt(mean)
A3 = rbind(c(T,T,T,T,T, T,T,T,F,F), c(T,T,T,T,T, F,F,F,F,F), c(T,T,F,F,F, F,F,F,F,F))
x3 = c(9,8,1,1,1, 1,1,1,1,1)
x3 = x3 / sum(x3)
b3 = A3 %*% x3
r3 = approx.region.gamma.1(A3, b3)

