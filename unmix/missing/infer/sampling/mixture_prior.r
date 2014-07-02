# Priors for the mixture proportions.

library(mnormt)

source("git/utils.r")

# A prior on the sort matrix.
# Args:
#   m - the sort matrix (with rows unnormalized); entries
#     should be between 0 and 1
#   m.cor - how correlated each sort matrix entry is with
#     the volume of each cell; this reflects how "noisy"
#     our estimate for the sorting of a particular cell is
#   vol.mean, vol.sd - prior on the volume of each cell
# Returns: a function which returns the log-likelihood
#   of the given sort matrix. The first row of the sort
#   matrix is presumed to be "all" (and so corresponds to
#   cell volume), while further rows correspond to the rows
#   of the sort matrix.
sort.prior.mvnorm.1 = function(m, m.cor, vol.mean, vol.sd) {
  # the dimension of the distribution for each cell
  n = 1 + nrow(m)

  # these will store the parameters of the normal distributions
  mu = list()
  si = list()

  for(j in 1:ncol(m)) {
    mu1 = rep(0, n)
    mu1[1] = vol.mean[
    si1 = matrix(0, nrow=n, ncol=n)
    si1[1,1] = vol.var[j]

    for(i in 1:nrow(m)) {
      mu
      si1[i+1, i+1] = vol.var[j]
      si1[1, i+1] = s

    }

  }

  list(mu = mu, si = si)
}

# Given the above parameters, computes the log-likelihood.
# Args:
#   p - parameters, including "mu" and "si" (e.g., from
#     sort.prior.mvnorm.1)
#   x - a sort matrix
# Returns: log-likelihood of that sort matrix.
sort.prior.ll = function(p) function(x) {
  function(x) {
    # Since we assume each cell is independent, we can sum up
    # the log-likelihood for each independently.
    s = 0
    for(j in 1:ncol(x)) {
      # XXX this is inefficient, as it keeps re-solving si
      s = s + dmnorm(x[,j], mean = p$mu[[j]], varcov=p$si[[j]],
        log=TRUE)
    }
    s
  }
}







