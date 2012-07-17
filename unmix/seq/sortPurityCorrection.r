# Correcting for sorting purity.

# get reads per million
r = as.matrix(read.table("git/unmix/seq/quant/readsPerMillion.tsv.gz",
  sep="\t", header=TRUE, check.names=FALSE, row.names=1))

# get sorting statistics
sorting.stats = read.table("git/unmix/seq/sorting.stats.tsv",
  sep="\t", header=TRUE, as.is=TRUE)
avg.sort.purity = c(by(sorting.stats$purity, sorting.stats$gene, mean))

# Does the correction for cases in which we have measured sort purity.
# Args:
#   r.mean, r.var - matrices of mean and variance of read counts
#   p - vector of sort purity
# Returns: r, with any fractions present in p corrected
correct.for.purity = function(r.mean, r.var, p) {
  for(g in colnames(r))
    if (g %in% names(p)) {
      r.mean[,g] = ( r.mean[,g] - (1-p[g]) * r.mean[,"all"] ) / p[g]
      r.var[,g] = ( r.var[,g] + (1-p[g]) * r.var[,"all"] ) / p[g]
    }

  r.mean[ r.mean < 0 ] = 0
  r.var[ r.var < 0 ] = 0

  list(r.mean = r.mean, r.var = r.var)
}

# correct this (note that we're assuming variance = mean)
r.corrected = correct.for.purity(r, r, avg.sort.purity)

