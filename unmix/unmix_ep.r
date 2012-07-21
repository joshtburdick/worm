# Unmixes using EP.






# Unmixes using EP.
# Args:
#   m - the cell-sorting matrix
#   x, x.var - mean and variance of the expression in each fraction,
#     as matrices with one row per gene, one column per fraction
#   output.dir - directory in which to write output
# Side effects: writes one .Rdata file per gene, each containing:
#   m, v - mean and variance of the posterior
#   x, x.var - the per-fraction mean and variance that were used
#   terms - the final term approximations
#   update.stats - how much the term approximations changed at each step
unmix.ep = function(m, x, x.var, output.dir) {






}

