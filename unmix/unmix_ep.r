# Unmixes using EP.

load("git/unmix/image/sort_matrix.Rdata")

# the read depth, with and without correction for sorting purity
source("git/unmix/seq/sortPurityCorrection.r")

source("git/unmix/ept/approx_region.r")

# scale rows of this to add up to 1
m = sort.matrix / apply(sort.matrix, 1, sum)

# limit to cases in which we have measurements
m = m[ colnames(r.corrected$r.mean) , ]


# Unmixes using EP.
# Args:
#   m - the cell-sorting matrix
#   x, x.var - mean and variance of the expression in each fraction,
#     as matrices with one row per gene, one column per fraction
#   output.dir - directory in which to write output
# Side effects: writes one .Rdata file per gene, each containing:
#   m, v - mean and variance of the posterior
#   x, x.var - the per-fraction mean and variance that were used
#   pos.terms - the final term approximations (which, with the prior, can
#     be used to reconstruct the full covariance posterior)
#   update.stats - how much the term approximations changed at each step
unmix.ep = function(m, x, x.var, output.dir) {
  system(paste("mkdir -p", output.dir))

  # restrict to cases in which the column sums to > 0
  nz = apply(m, 2, sum) > 0

  m1 = m[ , nz ]

  for(g in rownames(x)) {
    cat(g, "\n")
    r = approx.region(m1, x[g,], x.var[g,], prior.var=1e3)
    ep = list(m = r$m, v = r$v, x = x, x.var = x.var,
      pos.terms = r$pos.terms, update.stats = r$update.stats)
    save(ep, file = paste(output.dir, "/", g, ".Rdata", sep=""))
  }
}

# Summarizes all the results in one directory.
# Side effects:
summarize = function(result.dir) {


}

# a test
if (TRUE) {
  g = c("pha-4", "ceh-6", "irx-1")
  unmix.ep(m, r.corrected$r.mean[g,], r.corrected$r.var[g,],
    "git/unmix/ep.20120723")
}

