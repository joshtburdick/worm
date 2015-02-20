# Some tests of NMF.

source("git/unmix/ml/nmf.r")

# Slavish, slow translation of (5), for testing.
nmf.div.update.naive = function(W, H, V) {
  H1 = H
  WH = W %*% H

  z = apply(W, 2, sum)

  for(a in 1:nrow(H))
    for(m in 1:ncol(H)) {
      s = 0
      for(i in 1:nrow(W))
        s = s + W[i,a] * V[i,m] / WH[i,m]

      H1[a,m] = H[a,m] * s / z[a]
    }

  H1
}

nmf.test.1 = function(m, n, k) {
  W = matrix(1, nrow=m, ncol=k)
  H = matrix(1, nrow=k, ncol=n)
  V = matrix(rgamma(m*n, shape=0.5, rate=1), nrow=m, ncol=n)

  r = nmf.1(W, H, V, nmf.euclidean.update)
  r
}

nmf.test.2 = function(m, n, k) {
  W = matrix(rgamma(m*k, shape=1, rate=1), nrow=m, ncol=k)
  H = matrix(rgamma(k*n, shape=1, rate=1), nrow=k, ncol=n)
  V = W %*% H
  W = matrix(rgamma(m*k, shape=1, rate=1), nrow=m, ncol=k)
  H = matrix(rgamma(k*n, shape=1, rate=1), nrow=k, ncol=n)
  r = nmf.1(W, H, V, nmf.div.update, iters=100)
  r
}


r1 = nmf.test.1(3,5,7)

r2 = nmf.test.1(3,7,5)

r3 = nmf.test.2(3,5,7)
r4 = nmf.test.2(3,7,5)

