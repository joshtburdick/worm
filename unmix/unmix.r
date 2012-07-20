# Unmixing using the pseudoinverse.

load("git/unmix/image/sort_matrix.Rdata")

# the read depth, with and without correction for sorting purity
source("git/unmix/seq/sortPurityCorrection.r")

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
    A1 = M/sqrt(as.vector(b.var[g,]))
    B1 = b[g,]/as.vector(b.var[g,])

#    r = lsei(A = rbind(A1, 1e-6 * diag(ncol(M))), B = c(B1, rep(0, ncol(M))),
#      G = Diagonal(ncol(M)), H = rep(0, ncol(M)))
    r = lsei(A = A1, B = B1,
      G = Diagonal(ncol(M)), H = rep(0, ncol(M)))

    x[g,] = r$X
  }

  x
}

# scale rows of this to add up to 1
m = sort.matrix / apply(sort.matrix, 1, sum)

# limit to cases in which we have measurements
m = m[ colnames(r.corrected$r.mean) , ]

x.pseudoinverse = { 
  x = r.corrected$r.mean %*% pseudoinverse(t(m))
  x[,"P0"] = 0
  x
}


test1 = function() {
  x11()
  par(mfrow=c(5,1))
  avg.expr = apply(r1, 1, mean)
  set.seed(0)
#  genes = sample(names(avg.expr)[avg.expr > 100], 5)
  genes = c("pha-4", "ceh-26", "pal-1", "rgs-3")
  r = unmix.lsei(m1, r1[genes,], r1[genes,])
  for(g in genes) {
    plot(r[g,], type="h", main=g)
  }
  r
}

# foo = test1()
