# Unmixes using lsei() function.

source("git/unmix/unmix.r")
source("git/unmix/ept/approx_region.r")

output.dir = "git/unmix/"

# Unmixes using the constraints that x >= 0, and _exact_ constraint
# that Mx = b (b.var is ignored.)
# Args:
#   M - the cell-sorting matrix
#   b, b.var - the mean and variance of the expression in each fraction,
#     as matrices with one row per gene, and one column per fraction
# Returns: the estimated expression in each cell
unmix.lsei.eq = function(M, b, b.var, prior.var = 1e2) {
  source("git/unmix/ept/approx_region.r")

  x = matrix(0, nrow = nrow(b), ncol=ncol(M))
  rownames(x) = rownames(b)
  colnames(x) = colnames(M)

  for(g in rownames(b)) {
    cat(g, "")
    try({
      r = lsei(E = M, F = b[g,],
        G = diag(ncol(M)), H = rep(0, ncol(M)), type=1)
      x[ g, ] = r$X
    })
  }

  x
}

# Unmixes using the constraints that x >= 0.
# Args:
#   M - the cell-sorting matrix
#   b, b.var - the mean and variance of the expression in each fraction,
#     as matrices with one row per gene, and one column per fraction
#   prior.var - variance of prior, as a constant times the maximum variance
#     in any fraction (without this, result is spiky)
# Returns: the estimated expression in each cell
unmix.lsei = function(M, b, b.var, prior.var = 1e2) {
  source("git/unmix/ept/approx_region.r")

  x = matrix(0, nrow = nrow(b), ncol=ncol(M))
  rownames(x) = rownames(b)
  colnames(x) = colnames(M)

  for(g in rownames(b)) {
    cat(g, "")
    try({

      # compute prior, with the fuzzy constraint observed
      p = prior.var * max(b.var[g,], na.rm=TRUE)
      mv = mvnorm.2(rep(0, ncol(M)), rep(p, ncol(M)),
        M, b[g,], b.var[g,])
      C = chol( chol2inv(chol(mv$V)) )
      
      # estimate just the non-zero cells
      r = lsei(A = C, B = C %*% mv$m,
        G = diag(ncol(M)), H = rep(0, ncol(M)), type=2)

      # only write the cells in non-zero fractions
      x[ g, ] = r$X
    })
  }

  x
}

compute.lsei = function() {

  g = gene.list

  x.lsei = unmix.lsei(M, r.corrected$m[g,], r.corrected$v[g,], prior.var=1e2)
  save(x.lsei, file=paste(output.dir, "/x.lsei.Rdata", sep=""))

  x.lsei.time = unmix.lsei(M.facs.and.ts,
      cbind(r.corrected$m[g,], r.ts$m[g,]),
      cbind(r.corrected$v[g,], r.ts$v[g,]))
  save(x.lsei.time, file=paste(output.dir, "/x.lsei.time.Rdata", sep=""))
}

