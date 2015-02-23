# Trying using Lee & Seung's NMF method from
# "Algorithms for Non-negative Matrix Factorization".
# FIXME find a citation

source("git/utils.r")

# The Euclidean rule, from (4) in that paper.
# I'm only implementing the case for updating H, since
# updating W is symmetric.
# Args:
#   W, H, V - the matrices for which we want WH = V
# Returns: updated H
nmf.euclidean.update = function(W, H, V)
  H * ( (t(W) %*% V) / (t(W) %*% W %*% H) )

# The divergence update rule, from (5) in that paper.
# Again, only implementing one case.
# ??? do I need both cases?
# Args:
#   W, H, V - the matrices for which we want WH = V
# Returns: updated H
nmf.div.update = function(W, H, V) {
  A = t(W) %*% (V / (W %*% H))
  M = A / apply(W, 2, sum)
  H * M
}

# Does NMF.
# Args:
#   W, H, V - the matrices in question
#   update.rule - the update rule to use
#   iters - number of iterations to do
# Returns: list with components
#   W, H - NMF results
#   update.stats - how much each thing changed
nmf.1 = function(W, H, V, update.rule, iters=50) {
  update.stats = NULL

  for(iter in 1:iters) {
    H1 = update.rule(W, H, V)
    W1 = t( update.rule(t(H1), t(W), t(V)) )
    update.stats = rbind(update.stats,
      c(W = max(abs(W-W1)), H = max(abs(H-H1)),
        abs.diff = max(abs(W %*% H - V))))
    H = H1
    W = W1
  }

  list(W = W, H = H, update.stats = update.stats)
}

# Unmixing, mostly based on this NMF method.
# (This also enforces that certain things add up to
# 1, though I don't know if that's needed.)
# Args:
#   m - the sort matrix (with one row per sort fraction,
#     and one column per cell)
#   r - the read data (with one row per gene)
#   x - starting point for expression (optional)
#   max.iters - max. updates to do
#   eps - stop when updates are smaller than this
# Returns: a list with components
#   m - the modified sort matrix
#   x - the unmixed data
#   update.stats - stats about convergence
unmix.nmf = function(m, r, x=NULL, max.iters=100, eps=1e-10) {
  # XXX transposing this, because it's what I'm used to
  r = t(r)
  update.rule = nmf.div.update
  update.stats = NULL

  if (is.null(x)) {
    x = matrix(1, nrow=ncol(m), ncol=ncol(r))
    rownames(x) = colnames(m)
    colnames(x) = colnames(r)
  }

  for(iter in 1:max.iters) {
# browser()
    x1 = update.rule(m, x, r)
    x1 = x1 / apply(x1, 1, sum)
    m1 = t( update.rule(t(x1), t(m), t(r)) )
    m1 = m1 / apply(m1, 1, sum)

    s1 = c(x = max(abs(x-x1)), m = max(abs(m-m1)),
        abs.diff = max(abs(m %*% x - r)))
write.status(paste(iter, s1[1], s1[2], s1[3]))
    if (max(s1, na.rm=TRUE) <= eps)
      break
    update.stats = rbind(update.stats, s1)
    
    x = x1
    m = m1
  }

  rownames(update.stats) = 1:nrow(update.stats)
  list(m = m, x = t(x), update.stats = update.stats)
}

