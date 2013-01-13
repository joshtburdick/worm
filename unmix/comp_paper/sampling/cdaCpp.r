# Attempt at implementing the "cda" variant of the
# hit-and-run method.

library("Rcpp")

sourceCpp(file = "git/unmix/comp_paper/sampling/cda.cpp")


# Tiny test; samples three numbers adding up to one.
cdaCpp.test0 = function() {
  A = matrix(c(1,1,1)/3, ncol=3, byrow=TRUE)
  b = c(1)
#  x0 = c(1,1,1) / 3
  x0 = c(0.99,0.005,0.005)
  Z = matrix(c(2,-1,-1, -1,2,-1, -1,-1,2), ncol=3, byrow=TRUE)

  r = cdaCpp(A, b, x0, Z, 1000000, 0)

  png("git/unmix/comp_paper/sampling/cdaCppTest0.png",
    width=900, height=300)
  par(mfrow=c(1,3))
  for(j in c(1:3))
    hist(r[,j], breaks=100, col="grey")
  dev.off()
}


