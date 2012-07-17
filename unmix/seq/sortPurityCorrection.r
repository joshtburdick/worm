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
#   r - matrix of read counts
#   p - vector of sort purity
# Returns: r, with any fractions present in p corrected
correct.for.purity = function(r, p) {
  for(g in colnames(r))
    if (g %in% names(p))
      r[,g] = ( r[,g] - (1-p[g]) * r[,"all"] ) / p[g]

  r
}

r.corrected = correct.for.purity(r, avg.sort.purity)

