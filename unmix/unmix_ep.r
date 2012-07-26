# Unmixes using EP.

source("R/lineage/tree_utils.r")

load("git/unmix/image/sort_matrix.Rdata")

# the read depth, with and without correction for sorting purity
source("git/unmix/seq/sortPurityCorrection.r")

source("git/unmix/ept/approx_region.r")
# source("R/unmix/ep/ep.diag.3.r")

# scale rows of this to add up to 1
M = sort.matrix / apply(sort.matrix, 1, sum)

# limit to cases in which we have measurements
M = M[ colnames(r.corrected$r.mean) , ]


# Unmixes using EP.
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
    r = approx.region(M1, x[g,], x.var[g,], prior.var=Inf)    # 100 * max(x.var[g,]))
    ep = list(g = g, m = r$m, v = r$v, x = x[g,], x.var = x.var[g,],
      t = r$t, update.stats = r$update.stats)

#    ep = ep.diag.3(rep(1e-1, 1340), 1e3 * diag(1340),
#      m1, x[g,], b.var=x.var[g,],
#      max.iters=30, converge.tolerance=1e-5)

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

  num.cells = apply(cell.lineage.matrix, 1, sum)

  ep.summary = array(dim=c(length(genes), length(lin.node.names), 4),
    dimnames=c("gene", "cell", "stat"))
  dimnames(ep.summary)[[1]] = genes
  dimnames(ep.summary)[[2]] = colnames(M)
  dimnames(ep.summary)[[3]] =
    c("per.cell.mean", "per.cell.var", "lineage.mean.total", "lineage.var.total")

  for(i in 1:length(files)) {
    cat(genes[i], "")
    ep = NULL
    load(paste(result.dir, "/", files[i], sep=""))
    r = ep.get.result(ep, M)
    ep.summary[genes[i],,"per.cell.mean"] = r$m
    ep.summary[genes[i],,"per.cell.var"] = diag(r$V)
    ep.summary[genes[i],,"lineage.mean.total"] = cell.lineage.matrix %*% r$m / num.cells
    ep.summary[genes[i],,"lineage.var.total"] =
      diag(cell.lineage.matrix %*% r$V %*% t(cell.lineage.matrix)) / num.cells
  }

  save(ep.summary, file=output.file)
}

# these are for getting the names of genes to unmix
load("R/unmix/comp_paper/expr.cell.Rdata")
enriched.fraction = read.table("R/unmix/sort_paper/unmix/fraction/enriched.fraction.tsv",
  header=TRUE, sep="\t", row.names=1, as.is=TRUE)

# Runs EP on some genes
run.ep = function() {
  g = union(rownames(expr.cell), rownames(enriched.fraction))
  g = intersect(g, rownames(r.corrected$r.mean))
  g = g[1:20]
  unmix.ep(M, r.corrected$r.mean[g,], r.corrected$r.var[g,],
    "git/unmix/ep.20120726a")

#  ep.summarize.dir("git/unmix/ep.20120724/", "git/unmix/ep.20120724.summary.Rdata", M)
}


