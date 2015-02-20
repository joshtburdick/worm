# Trying using Lee & Seung's NMF method from
# "Algorithms for Non-negative Matrix Factorization".
# FIXME find a citation

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

  for(i in 1:iters) {
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






