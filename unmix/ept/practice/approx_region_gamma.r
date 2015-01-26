# Another attempt at approximating marginals with gamma
# distributions.

source("git/unmix/ept/gamma.r")


# wrapper for this
gca = function(x, A, b) {
  x1 = array(x, dim=c(dim(x)[1], dim(x)[2], 1))
  dimnames(x1)[[1]] = c("e1", "e2")
  x2 = gamma.conditional.approx(x1, A, b)
  x2[,,1]
}


# Approximates marginals with a gamma distribution.
# Not focussing on efficiency here...
# Args:
#   A, b - these give the constraints
#   num.iters - number of iterations
# Returns: list with elements
#   x - the posterior
#   term - the term approximations
#   FIXME update stats?
# XXX not working
approx.region.gamma = function(A, b, num.iters=40) {

  # approximating terms
  term = array(1, dim=c(2, ncol(A), nrow(A)))
  dimnames(term)[[1]] = c("m", "v")
  term = gamma.mv2n(term)

  # prior
  prior = 1 * gamma.mv2n(rbind(m=rep(1,ncol(A)), v=rep(1,ncol(A))))

  # the posterior  
  x = prior + apply(term, c(1,2), sum)

  # iterate
  for(iter in 1:num.iters) {

    # compute message to each term
    t1 = term
    for(i in 1:nrow(A)) {
      t1[,,i] = x - term[,,i]
    }

    # moment-match
    x1 = gamma.conditional.approx(t1, A, b)

    # correct terms
    for(i in 1:nrow(A)) {
      term[,,i] = term[,,i] + 1 * (x1[,,i] - x)
    }

    # update posterior    
    x = prior + apply(term, c(1,2), sum)
  }
  list(x = x, term = term)
}

# test of how this converges
t1 = function() {
  A = matrix(rgamma(1000, shape=1, scale=1), nrow=1)
  b = 1
  t0 = gamma.mv2n(rbind(m=rep(1,1000), v=rep(1,1000)))

  q = gca(t0, A, b);  gamma.n2mv(q[,1:5]); t0 = q - t0
  print(q)

}

# more trying to get a simple case of this working
t2 = function() {
  f = function(x) gca(x, t(c(1,2,4)), 1)

#  q0 = gamma.mv2n(rbind(m=c(1,1,1), v=c(1,1,1)))
  t0 = gamma.mv2n(rbind(m=c(1,1,1), v=c(1,1,1)))
  q = f(t0)
  for(i in 1:20) {
    t0 = q - t0
    q = f(t0)
    print(gamma.n2mv(q))
  }
}



# a small test
A2 = matrix(c(1,1,1,1,0,0, 0,0,0,1,1,1), nrow=2, byrow=TRUE)
A2 = rbind(A2, 1-A2)
# b2 = c(0.6, 0.8)
# b2 = c(b2, 1-b2)
# x2 = array(gamma.s2n(rbind(a=rep(1,6), b=rep(1,6))), dim=c(2,6,1))
# dimnames(x2)[[1]] = c("e1", "e2")

x2 = rgamma(6, 2, 1)
x2 = x2 / sum(x2)
b2 = as.vector(A2 %*% x2)


foo = approx.region.gamma(A2, b2)





