# Unmixes using EP.

source("git/unmix/unmix.r")
source("git/unmix/ept/approx_region.r")

# Unmixes some genes using EP.
# Args:
#   M - the cell-sorting matrix
#   x, x.var - mean and variance of the expression in each fraction,
#     as matrices with one row per gene, one column per fraction
#   output.dir - directory in which to write output
# Side effects: writes one .Rdata file per gene, each containing:
#   m, v - mean and variance of the posterior
#   x, x.var - the per-fraction mean and variance that were used
#   pos.terms - the final term approximations (which, with the prior, can
#     be used to reconstruct the full covariance posterior)
#   update.stats - how much the term approximations changed at each step
unmix.ep = function(M, x, x.var, output.dir) {
  system(paste("mkdir -p", output.dir))

  # restrict to cases in which the column sums to > 0
  nz = apply(M, 2, sum) > 0
  M1 = M[ , nz ]

  for(g in rownames(x)) {
    cat("\n", g, "")
    try({
      r = approx.region(M, x[g,], x.var[g,], prior.var=1e3 * max(x.var[g,]))
      ep = list(g = g, m = r$m, v = r$v, x = x[g,], x.var = x.var[g,],
        t = r$t, update.stats = r$update.stats)
      save(ep, file = paste(output.dir, "/", g, ".Rdata", sep=""))
    })
  }
}

# Given an EP result:
# - expands the full covariance matrix (by adding in the "Ax ~ b" constraint)
# - adds in columns that were removed (because of zeros in cell-sorting matrix)
# Args:
#   ep - an EP result, including "pos.terms", "x", and "x.var"
#   M - the cell-sorting matrix that was used
# Returns: m and V, the mean and covariance of the posterior
ep.get.result = function(ep, M) {

  # restrict to cases in which the column sums to > 0
  nz = apply(M, 2, sum) > 0
  M1 = M[ , nz ]

  # get full covariance matrix
  r = mvnorm.2(ep$t[,"m"], ep$t[,"v"], M1, ep$x, ep$x.var) 

  # expand to include things that were zero
  m = rep(0, ncol(M))
  names(m) = colnames(M)
  V = matrix(0, nrow=ncol(M), ncol=ncol(M))
  rownames(V) = colnames(M)
  m[nz] = r$m
  V[nz,nz] = r$V

# correct for cell volume?
#  m = m / as.vector(M["all",])
#  V = V / as.vector(M["all",])
#  V = t( t(V) / as.vector(M["all",] ))

  # set NA cases to 0
  m[ is.na(m) ] = 0
  V[ is.na(V) ] = 0

  list(m = m, V = V)
}

# Summarizes all the results in one directory.
# Args:
#   result.dir - directory of results to summarize
#   output.file - where to write output to
#   M - the original sorting matrix which was used
# Side effects: writes an .Rdata file containing an array "ep.summary",
#   with dimensions "gene", "cell", and "stat", where "stat" includes
#   mean and variance of the per-cell and lineage totals
ep.summarize.dir = function(result.dir, output.file, M) {
  files = sort(list.files(result.dir, pattern=".Rdata"))
  genes = sub(".Rdata$", "", files)

  # these are restricted to just the cells which were unmixed
  num.cells = apply(cell.lineage.matrix, 1, sum)[colnames(M)]
  clm = cell.lineage.matrix[ colnames(M), colnames(M) ]

  ep.summary = array(0, dim=c(length(genes), length(lin.node.names), 4),
    dimnames=c("gene", "cell", "stat"))
  dimnames(ep.summary)[[1]] = genes
  dimnames(ep.summary)[[2]] = lin.node.names
  dimnames(ep.summary)[[3]] =
    c("per.cell.mean", "per.cell.var", "lineage.mean.total", "lineage.var.total")

  for(i in 1:length(files)) {
    cat(genes[i], "")
    ep = NULL
    load(paste(result.dir, "/", files[i], sep=""))
    r = ep.get.result(ep, M)
    ep.summary[genes[i],colnames(M),"per.cell.mean"] = r$m
    ep.summary[genes[i],colnames(M),"per.cell.var"] = diag(r$V)
    ep.summary[genes[i],colnames(M),"lineage.mean.total"] = clm %*% r$m / num.cells
    ep.summary[genes[i],colnames(M),"lineage.var.total"] =
      diag(clm %*% r$V %*% t(clm)) / num.cells
  }

  save(ep.summary, file=output.file)
}

run.ep = function() {
  g = gene.list
  unmix.ep(M, r.corrected$m[g,], r.corrected$v[g,],
    "git/unmix/ep.facs.only")
  unmix.ep(M.facs.and.ts,
    cbind(r.corrected$m[g,], r.ts$m[g,]),
    cbind(r.corrected$v[g,], r.ts$v[g,]),
    "git/unmix/ep.with.time")

  ep.summarize.dir("git/unmix/ep.facs.only/", "git/unmix/ep.facs.only.Rdata", M)
  ep.summarize.dir("git/unmix/ep.with.time/", "git/unmix/ep.with.time.Rdata", M.facs.and.ts)
}

