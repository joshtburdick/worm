# Another attempt at approximating marginals with gamma
# distributions.

source("git/unmix/ept/gamma.r")
source("git/unmix/ept/practice/sd.r")
source("git/unmix/ept/practice/gamma_conditional_6.r")
source("git/unmix/ept/practice/gamma_conditional_numerical.r")

# wrapper for this
gca = function(x, A, b) {
  x1 = array(x, dim=c(dim(x)[1], dim(x)[2], 1))
  dimnames(x1)[[1]] = c("e1", "e2")
  x2 = gamma.conditional.approx(x1, A, b)
  x2[,,1]
}

# Approximates marginals of
#   x_i ~ Gamma() | ax = b
# Args:
#   x - gamma distribution (natural parameters)
#   a, b - the linear constraint (for now, assuming
#     all a != 0)
# Returns: natural parameters of x | Ax = b
# FIXME: compare this with sampling?
# XXX deprecated; see gamma_conditional_numerical.r .
gca1 = function(x, a, b) {
  # convert to moments
  m = gamma.n2mv(x)

  # scale this
  s = b / a
  m["m",] = m["m",] / a
  m["v",] = m["v",] / (a*a)

  # approximate marginals, if those sum to 1
  #  r = sd.s2mv(gamma.mv2s(m))
  #  r = gamma.n2mv(gamma.cond.sum1(gamma.mv2n(m)))
  r = gamma.cond.sum.numerical.int.1(gamma.mv2n(m))

  # undo scaling
  r["m",] = r["m",] * s
  r["v",] = r["v",] * (s*s)

  # convert back to natural parameters
  gamma.mv2n(r)
}

# Approximates marginals with a gamma distribution.
# Not focussing on efficiency here...
# Args:
#   A, b - these give the constraints
#   num.iters - number of iterations
# Returns: list with elements
#   x - the posterior
#   term - the term approximations
#   x.log - list of posterior approximations
#   FIXME update stats?
# XXX not sure how well this is working
approx.region.gamma = function(A, b, num.iters=50) {
  x.log = list()

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
    for(i in 1:nrow(A)) {
      t1 = x - term[,,i]
#      x1 = gamma.cond.sampling(t1, t(as.vector(A[i,])), b[i])
      x1 = gamma.cond.sum.numerical(t1, A[i,], b[i])
#      x1 = gamma.cond.orig(t(A[i,]), b[i])(t1)
      term[,,i] = 1 * term[,,i] + 0.2 * (x1 - x)
#      term[,,i] = x1 - t1
    }

    # update posterior
    x = prior + apply(term, c(1,2), sum)
    x.log[[iter]] = x
print(gamma.n2mv(x))
  }

  list(x = x, term = term, x.log = x.log)
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
t3 = function() {

  A2 = matrix(c(1,1,0.1, 0.1,1,1), nrow=2, byrow=TRUE)
  # A2 = rbind(A2, 1-A2)
  # b2 = c(0.6, 0.8)
  # b2 = c(b2, 1-b2)
  # x2 = array(gamma.s2n(rbind(a=rep(1,6), b=rep(1,6))), dim=c(2,6,1))
  # dimnames(x2)[[1]] = c("e1", "e2")

  x2 = c(1,1,1)
  b2 = as.vector(A2 %*% x2)
  approx.region.gamma(A2, b2)
}

t4 = function() {
  A = matrix(sample(0:4, 400, replace=TRUE)+0.1, nrow=5)
  x = rgamma(80, 0.5, 10)
  b = A %*% x
  r = approx.region.gamma(A, b)
  list(A = A, x = x, b = b, r = r, y = gamma.n2mv(r$x)["m",])
}


