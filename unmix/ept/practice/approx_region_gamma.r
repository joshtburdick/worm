# Another attempt at approximating marginals with gamma
# distributions.

source("git/unmix/ept/gamma.r")


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
  prior = 1 * gamma.mv2n(rbind(m=rep(10,ncol(A)), v=rep(100,ncol(A))))

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
      term[,,i] = term[,,i] + x1[,,i] - x
    }

    # update posterior    
    x = prior + apply(term, c(1,2), sum)
  }

  list(x = x, term = term)
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







