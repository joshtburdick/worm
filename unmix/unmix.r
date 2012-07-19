# Unmixing using the constraints.




# Unmixes using the constraints.
# Args:
#   M - the cell-sorting matrix
#   b, b.var - the mean and variance of the expression in each fraction,
#     as matrices with one row per gene, and one column per fraction
# Returns: the estimated expression in each cell
unmix.lsei = function(M, b, b.var) {
  x = matrix(nrow = nrow(b), ncol=ncol(M))
  rownames(x) = rownames(b)
  colnames(x) = colnames(M)

  for(g in rownames(b)) {
    cat(g, "")

# XXX this gives a very spiky answer
#    r = lsei(A = M/sqrt(as.vector(b.var[g,])), B = b[g,]/as.vector(b.var[g,]),
#      G = Diagonal(ncol(M)), H = rep(0, ncol(M)))

# ... and this is saying the inequalities are contradictory
    r = lsei(E = M / sqrt(as.vector(b.var[g,])), F = b[g,] / as.vector(b.var[g,]),
      G = Diagonal(ncol(M)), H = rep(0, ncol(M)))
    x[g,] = r$X
  }

  x
}




