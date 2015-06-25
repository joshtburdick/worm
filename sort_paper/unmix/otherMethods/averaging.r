# Trying to use a weighted average of per-cell measurements.

source("git/sort_paper/unmix/otherMethods/crossValAccuracyRpm.r")

# Unmixes using a weighted average of all the measurements
# for a given cell.
unmix.average1 = function(m, r) {
  m.columns.normalized = t( t(m) / apply(m, 2, sum) )
  x = r %*% m.columns.normalized
  list(x = x)
}

write.crossval.graphs(unmix.average1, "average1")

