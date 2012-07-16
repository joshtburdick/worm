# Factor which forces a variable to be positive.







# A factor which forces a variable to be positive.
# Args:
#   x - the variable in question
positive.factor = function(x) {
  a = list
  a$update = function(m) {

  }
  # this is the log of the proportion of the distribution that's positive
  a$log.evidence = function(m) {
    b = canonical.to.standard.params(m$b)
    sum(pnorm(0, -b$mu, sd, log.p=TRUE))
  }
  a
}


