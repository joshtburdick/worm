# Does unmixing using the pseudoinverse.

library(corpcor)

source("unmix_test.r")

# Does unmixing of one gene.
# Function which, given a sort matrix and expression
#   in each fraction, returns unmixed expression.
unmix.pseudoinverse = function(m, x.fraction) {

  # XXX recomputing the pseudoinverse each time,
  # which is somewhat wasteful
  x = pseudoinverse(m) %*% as.vector(x.fraction)
  x = as.vector( x[,1] )
  x[ x < 0 ] = 0
  list(x = x)
}

run.unmix(unmix.pseudoinverse, "pseudoinverse/")

