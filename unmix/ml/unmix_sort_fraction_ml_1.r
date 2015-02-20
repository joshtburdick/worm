# Another attempt at estimating expression, and the sort matrix,
# simultaneously.

source("git/utils.r")
source("git/unmix/ml/pos_linear_solve.r")

# Estimates expression, and the sort matrix.
# Args:
#   m - sort matrix estimate (each row should sum to 1)
#   r - read data, as proportions (each column should sum to 1)
#   max.iters - maximum number of iterations to do
# Returns: list with elements:
#   m - modified sort matrix
#   x - mean of unmixed expression
#   update.stats - how much m and x changed
unmix.expr.and.sort.matrix.1 = function(m, r, max.iters=20) {
  update.stats = NULL
#  x = pos.linear.solve(m, r, max.iters=5, normalize="rows")$X
  x = 0

  for(iter in 1:max.iters) {
cat(paste0("\niter = ", iter, "\n"))
    x1 = pos.linear.solve(m, r, max.iters=5, eps=1e-10, normalize="rows")$X
# browser()
    m1 = t( pos.linear.solve(t(x1), t(r), max.iters=5, eps=1e-10, normalize="columns")$X )

    update.stats.1 = c(x = max(abs(x1-x)), m = max(abs(m1-m)))
write.status(paste(update.stats.1, collapse=" "))
    update.stats = rbind(update.stats, update.stats.1)
    m = m1
    x = x1
  }

  list(m = m, x = x, update.stats = update.stats)
}

