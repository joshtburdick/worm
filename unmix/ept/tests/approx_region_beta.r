# some tests



source("git/unmix/ept/approx_region_beta.r")



# slowish version of lin.constraint, which loops through
# each constraint separately.
lin.constraint.slow = function(m, v, A, b) {
  n = nrow(m)
  r.m = NULL
  r.v = NULL

  for(i in 1:n) {
    r = lin.constraint.1(m[i,], v[i,], A[i,,drop=FALSE], b[i], 0)
    r.m = rbind(r.m, r[,"m"])
    r.v = rbind(r.v, r[,"v"])
  }

  list(m = r.m, v = r.v)
}


compare.1 = function(m, v, A, b) {

  r1 = lin.constraint.slow(m, v, A, b)
#  print(r1)

  r2 = lin.constraint.1eq(m, v, A, b)
#  print(r2)

  cat("max. diff m =", max(abs(r1$m-r2$m)),
    "   v =", max(abs(r1$v-r2$v)), "\n")
}



# compare.1(c(3,4,1), c(1,2,3), t(c(1,1,1)), 1)

# compare.1(runif(10), runif(10), t(rep(1, 10)), 1)


m = matrix(1:6, nrow=2)
v = matrix(4:9, nrow=2)
A = matrix(1:6 %% 4, nrow=2)
b = c(2,3)

# print(lin.constraint.slow(m, v, A, b))
compare.1(m, v, A, b)

p1 = list(m = matrix(runif(15), nrow=3), v = matrix(runif(15), nrow=3),
  A = matrix(runif(15), nrow=3), b=c(7:9))
compare.1(p1$m, p1$v, p1$A, p1$b)

