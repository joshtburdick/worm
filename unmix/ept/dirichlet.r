# The Dirichlet distribution.

# Convert "usual" parameters to natural parameters.
dirichlet.p2n = function(b) b - 1

# Convert natural parameters to the "usual".
dirichlet.n2p = function(a) a + 1

# The (log) normalization constant, according to Wikipedia.
# (a.k.a. the "log-partition function").
# Args:
#   a - the natural parameters
logPartition.dirichlet = function(a)
  sum(lgamma(a+1)) - lgamma(sum(a+1))

# Log-evidence contribution for a Dirichlet variable.
# Args: list of Dirichlet distributions (in canonical parameters.)
log.evidence.dirichlet = function(a) {
  n = length(a)

  # parameters for the product of all of these
  s = 0 * a[[1]]
  for(i in 1:n)
    s = s + a[[i]]
  s = s

  dirichlet.log.partition(s) - sum(sapply(a, dirichlet.log.partition))
}


# A Dirichlet equality factor. ??? is this correct?
dirichlet.eq.factor = function(a) {
  x1 = a[[1]]
  x2 = a[[2]]

  list(update = list(x2, x1),
    log.evidence = log.partition.dirichlet(list(x1,x2)))
}

