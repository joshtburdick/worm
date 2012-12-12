
source("git/unmix/ept/gamma.r")

# Estimates marginals of gamma-distributed variables,
# conditional on a weighted sum having a gamma distribution.
# Args:
#   a - parameters of gamma variables
#   m - weighting vector
#   b - parameters of observation
# Returns: samples, with weights. 
cond.gamma.sample = function(a, m, b, num.samples=1000000) {
  x = NULL
  for(i in 1:ncol(a))
    x = cbind(x, rgamma(num.samples, shape=a["a",i], rate=a["b",i]))

# print(dim(x))
  s = t( m * t(x) )

  w = dgamma(s, shape=b["a",1], rate=b["b",1])
  x1 = apply(w*x, 2, sum) / sum(w)
  x2 = apply(w*x*x, 2, sum) / sum(w)
# print(length(x1))
# print(length(x2))
  si2 = x2 - x1^2
#  si2 = apply(w * t( t(x) - x1), 2, sum) / sum(w)

  list(p = cbind(m = x1, v = si2), n = sum(w>=0.8))  #, x = x, w = w)
}


a = gamma.mv2s(rbind(m=c(1,1,1), v=c(2,1,1)))
b = gamma.mv2s(rbind(m=1, v=1e-6))
# r = cond.gamma.sample(a, c(1,1,1), b)

# Does a bunch of sampling of these.
cond.gamma.sample.many = function() {
  cond.gamma.samples = NULL

  for(iter in 1:1000000) {
    a = matrix(rgamma(6, 1, 1), nrow=2, ncol=3)
    rownames(a) = c("a", "b")
    r = cond.gamma.sample(a, c(1,1,1), b, num.samples=1000000)
    r[["a"]] = a
    r[["p"]] = t(r[["p"]])

    cond.gamma.samples[[iter]] = r
    if (iter %% 10 == 0) {
      cat(iter, "")
      save(cond.gamma.samples,
        file="git/unmix/ept/practice/gamma_conditional_samples.Rdata")
    }
  }
}


