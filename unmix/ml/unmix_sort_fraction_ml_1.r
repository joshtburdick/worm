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
#   a.prior - prior on the sort matrix
#   x.prior - prior on expression
#   b - read.data, as proportions (each row should sum to 1)
#   max.iters - maximum number of iterations to do
#   save.hist - whether to save a history of estimates
# Returns: list with elements:
#   a - estimated sort matrix
#   x - estimated expression
#   update.stats - likelihood stats, and how much a and x changed
unmix.expr.and.sort.matrix.1 =
    function(a.prior, x.prior, b, max.iters=50, save.hist=FALSE) {

  # these track convergence (or lack thereof)
  update.stats = NULL
  h = list()

  # initialize estimates to the mean
  a = a.prior$a / a.prior$b

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

