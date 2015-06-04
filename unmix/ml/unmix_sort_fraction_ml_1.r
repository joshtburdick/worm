# Another attempt at estimating expression, and the sort matrix,
# simultaneously.

source("git/utils.r")
source("git/unmix/ml/pos_linear_solve.r")

# Does one step of "updating".
# Args:
#   A, B - these give the linear constraint
#   X - this gives the current solution
# Returns: an updated X, which is closer to a solution of
#   AX = B, with X >= 0, and each row of X summing to 1.
pos.linear.solve.1 = function(A, B, X) {
  X1 = lin.constraint(A, B)(X)
  X = move.pos(X, X1) + 1e-20
  X = X / apply(X, 1, sum)
  X
}

# Estimates expression, and the sort matrix.
# Args:
#   m - sort matrix estimate (each row should sum to 1)
#   r - read data, as proportions (each column should sum to 1)
#   max.iters - maximum number of iterations to do
#   save.x.history - whether to return all estimates of x
# Returns: list with elements:
#   m - modified sort matrix
#   x - mean of unmixed expression
#   update.stats - how much m and x changed
unmix.expr.and.sort.matrix.1 = function(m, r, max.iters=50, save.x.history=FALSE) {
  update.stats = NULL
  x.history = list()
#  x = pos.linear.solve(m, r, max.iters=5, normalize="rows")$X
#  x = 0
  x0 = matrix(1, nrow=ncol(m), ncol=ncol(r))
  x = pos.linear.solve.1(m, r, x0)  

  for(iter in 1:max.iters) {
cat(paste0("\niter = ", iter, "\n"))
#    x1 = pos.linear.solve(m, r, max.iters=5, eps=1e-10, normalize="rows")$X
#    m1 = t( pos.linear.solve(t(x1), t(r), max.iters=5, eps=1e-10, normalize="columns")$X )
# ??? why was I normalizing columns? that seems incorrect

    # XXX trying strictly alternating between optimizing these
    m1 = t(pos.linear.solve.1(t(x), t(r), t(m)))
    x1 = pos.linear.solve.1(m1, r, x)

    update.stats.1 = c(x = max(abs(x1-x)), m = max(abs(m1-m)))
write.status(paste(update.stats.1, collapse=" "))
    update.stats = rbind(update.stats, update.stats.1)
#    m = m1
    x = x1
    if (save.x.history) {
      x.history[[ iter ]] = t(x)
    }
  }

  list(m = m, x = x, update.stats = update.stats, x.history=x.history)
}

