# Unmixes using EP.

load("git/unmix/image/sort_matrix.Rdata")

# the read depth, with and without correction for sorting purity
source("git/unmix/seq/sortPurityCorrection.r")

source("git/unmix/ept/approx_region.r")
# source("R/unmix/ep/ep.diag.3.r")

# scale rows of this to add up to 1
m = sort.matrix / apply(sort.matrix, 1, sum)

# limit to cases in which we have measurements
m = m[ colnames(r.corrected$r.mean) , ]


# Unmixes using EP.
# Args:
#   m - the cell-sorting matrix
#   x, x.var - mean and variance of the expression in each fraction,
#     as matrices with one row per gene, one column per fraction
#   output.dir - directory in which to write output
# Side effects: writes one .Rdata file per gene, each containing:
#   m, v - mean and variance of the posterior
#   x, x.var - the per-fraction mean and variance that were used
#   pos.terms - the final term approximations (which, with the prior, can
#     be used to reconstruct the full covariance posterior)
#   update.stats - how much the term approximations changed at each step
unmix.ep = function(m, x, x.var, output.dir) {
  system(paste("mkdir -p", output.dir))

  # restrict to cases in which the column sums to > 0
  nz = apply(m, 2, sum) > 0

  m1 = m[ , nz ]

  for(g in rownames(x)) {
    cat(g, "\n")
try({
    r = approx.region(m1, x[g,], x.var[g,], prior.var=100 * max(x.var[g,]))
    ep = list(m = r$m, v = r$v, x = x, x.var = x.var,
      pos.terms = r$pos.terms, update.stats = r$update.stats)

#    ep = ep.diag.3(rep(1e-1, 1340), 1e3 * diag(1340),
#      m1, x[g,], b.var=x.var[g,],
#      max.iters=30, converge.tolerance=1e-5)

    save(ep, file = paste(output.dir, "/", g, ".Rdata", sep=""))
})
  }
}

# Given an EP result:
# - expands the full covariance matrix (by adding in the "Ax ~ b" constraint)
# - scales by cell size, and
# - and adds columns that were removed (because of zeros in cell-sorting matrix)
ep.summarize = function(ep, m) {




  list(m =  , V = )
}

# Summarizes all the results in one directory.
# Args:
#   result.dir - directory of results to summarize
#   output.file - where to write output to
#   m - the original sorting matrix which was used
# Side effects: writes an .Rdata file containing an array "ep.summary",
#   with dimensions "gene", "cell", and "stat", where "stat" includes
#   mean and variance of the per-cell and lineage totals
ep.summarize.dir = function(result.dir, output.file, m) {
  files = sort(list.files(result.dir, pattern=".Rdata"))
  genes = sub(".Rdata$", "", files)

  ep.summary = array(dim=c(nrow(expr.cell), length(lin.node.names), 4),
    dimnames=c("movie", "cell", "stat"))
  dimnames(a)[[1]] = genes
  dimnames(a)[[2]] = colnames(m)
  dimnames(a)[[3]] =
    c("per.cell.mean", "per.cell.var", "lineage.mean.total", "lineage.var.total")






}


# these are for getting the names of genes to unmix
load("R/unmix/comp_paper/expr.cell.Rdata")
enriched.fraction = read.table("R/unmix/sort_paper/unmix/fraction/enriched.fraction.tsv",
  header=TRUE, sep="\t", row.names=1, as.is=TRUE)

# a test
if (TRUE) {
  g = union(rownames(expr.cell), rownames(enriched.fraction))
  g = intersect(g, rownames(r.corrected$r.mean))
  g = g[1:20]
#  g = c("cnd-1", "cwn-1")
  unmix.ep(m, r.corrected$r.mean[g,], r.corrected$r.var[g,],
    "git/unmix/ep.20120724b")
}

