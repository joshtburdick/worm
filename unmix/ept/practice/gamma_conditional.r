
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
  for(i in 1:nrow(a))
    x = cbind(x, rgamma(num.samples, shape=a[i,"a"], rate=a[i,"b"]))

# print(dim(x))
  s = t( m * t(x) )

  w = dgamma(s, shape=b[1,"a"], rate=b[1,"b"])
  x1 = apply(w*x, 2, sum) / sum(w)
  x2 = apply(w*x*x, 2, sum) / sum(w)
# print(length(x1))
# print(length(x2))
  si2 = x2 - x1^2
#  si2 = apply(w * t( t(x) - x1), 2, sum) / sum(w)

  list(p = cbind(m = x1, v = si2), n = sum(w>=0.8))  #, x = x, w = w)
}




a = cbind(gamma.mv2s(cbind(m=c(2,1), v=c(1,1))))
# a = cbind(a=c(1,2), b=c(1,1))
b = gamma.mv2s(cbind(m=1, v=1e-6))

r = cond.gamma.sample(a, c(2,1), b)


