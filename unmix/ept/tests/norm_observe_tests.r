
source("git/unmix/ept/norm_observe.r")

A = matrix(c(1:6), nrow=2)

# some simple tests, in order of increasing weirdness
r = signif(lin.constraint.1(rep(0, 3), rep(1, 3), A, rep(1,2), rep(1,2)), 3)
stopifnot(!any(is.na(r)))
r = signif(lin.constraint.1(rep(0, 3), rep(1, 3), A, rep(1,2), rep(0,2)), 3)
stopifnot(!any(is.na(r)))
r = signif(lin.constraint.1(rep(0, 3), rep(1, 3), A, rep(1,2), rep(-1,2)), 3)
stopifnot(!any(is.na(r)))

# this seems to be well-defined, even when the variance is very high
r = signif(lin.constraint.1(rep(0, 3), c(1,1,1e5), A, rep(1,2), rep(0,2)), 3)
stopifnot(!any(is.na(r)))
r = signif(lin.constraint.1(rep(0, 3), c(1,1,1e5), A, rep(1,2), rep(1,2)), 3)
stopifnot(!any(is.na(r)))

# simpler cases, which may or may not fit my head...
r = lin.constraint.1(c(0,0,0), c(1,1,1), t(c(1,1,1)), 1, 0)
round(r, 3)

r = lin.constraint.1(c(0,0,0), c(1,1,1), t(c(1,1,1)), 1, 1)
round(r, 3)



r = lin.constraint.1(c(0,0,0), c(1,1,1e20), t(c(1,1,1)), 1, 1)
round(r, 3)

# note that if I simply drop that coordinate, I get
# a variance which is too small (not to mention a bogus mean)
r = lin.constraint.1(c(0,0), c(1,1), t(c(1,1)), 1, 1)
round(r, 3)

r = lin.constraint.1(rep(0,10), c(rep(1,9),1e8), t(rep(1,10)), 1, 1)
round(r, 3)




A2 = matrix(runif(20), nrow=2)
r = lin.constraint.1(rep(0,10), c(rep(0,9),1e5), A2, c(1,2), c(2,1))
round(r, 3)
r = lin.constraint.1(rep(0,10), c(rep(2,4),1e5,1e5,1e5,rep(3,3)), A2, c(1,2), c(2,1))
round(r, 3)


A3 = matrix(runif(30), nrow=3)
r = lin.constraint.1(rep(0,10), c(rep(0,9),1e5), A3, c(1,2,3), c(3,2,1))
round(r, 3)

r = lin.constraint.1(c(1:10), c(rep(1,4),1e8, 1e8, rep(10,4)), A3,
  c(1,2,3), c(2,3,1))
round(r, 3)

# Empirically, it looks like this would converge so that the things
# with finite variance are unaffected by the constraint, and the entries
# with infinite variance get determined by the constraint (plus the
# things with finite variance.) or something.
lin.constraint.1(c(0,0,0), c(1,2,1e5), t(c(1,1,1)), 2, 10)
lin.constraint.1(c(0,0,0), c(1,2,1e5), t(c(1,1,1)), 2, 1)
lin.constraint.1(c(0,0,0), c(1,2,1e9), t(c(1,1,1)), 2, 1)
lin.constraint.1(c(0,0,0), c(1,2,1e9), t(c(1,1,1)), 2, 100)

# Again, empiricallly, it seems that if there are at least as many
# Inf variance entries as constraints, then the constraints
# don't affect the entries with finite variance. 
r = lin.constraint.1(c(10:1), c(1:8,rep(1e7,2)), A3, c(1,2,3), c(0,0,0))
r = lin.constraint.1(c(10:1), c(1:7,rep(1e7,3)), A3, c(1,2,3), c(0,0,0))
r = lin.constraint.1(c(10:1), c(1:6,rep(1e7,4)), A3, c(1,2,3), c(0,0,0))


r = lin.constraint.1(c(10:3), c(1:8), A3[,c(1:8)], c(1,2,3), c(0,0,0))




