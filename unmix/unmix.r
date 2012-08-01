# Unmixing using the pseudoinverse.

library(corpcor)
library(limSolve)

load("git/unmix/image/sort_matrix.Rdata")

# the read depth, with and without correction for sorting purity
source("git/unmix/seq/sortPurityCorrection.r")

# limit sort matrix to only include things up to a certain time
t.cutoff = 420
cell.time.on.off =
  read.csv(gzfile("data/image/cellTimeOnOff.csv.gz"),
    row.names=1, as.is=TRUE)

# which cells to include
# XXX note that this is not in lineage order
cells.to.include =
  rownames(cell.time.on.off)[ cell.time.on.off$on < t.cutoff ]
M = sort.matrix[ , cells.to.include ]

# scale rows of this to add up to 1
M = M / apply(M, 1, sum)

# limit to cases in which we have measurements
M = M[ colnames(r.corrected$r.mean) , ]

# Removes cells in any fractions having zero expression.
# Args:
#   M - the sort matrix
#   b - expression in each fraction
# Returns: a list with components
#   M, b, b.var - the corresponding variables, with zeros
#     removed from b (and M and b.var; M will also have
#     columns corresponding to zero fractions removed)
#   nz - which columns were zero
remove.zeros = function(M, b, b.var) {
  z.fraction = (b == 0)
  z = apply(M[z.fraction,], 2, sum) > 0
  list(M = M[ !z.fraction, !z ], b = b[!z.fraction], b.var = b.var[!z.fraction],
    nz = !z)
}

# Unmixes using the constraints that x >= 0.
# Args:
#   M - the cell-sorting matrix
#   b, b.var - the mean and variance of the expression in each fraction,
#     as matrices with one row per gene, and one column per fraction
# Returns: the estimated expression in each cell
unmix.lsei = function(M, b, b.var) {
  x = matrix(0, nrow = nrow(b), ncol=ncol(M))
  rownames(x) = rownames(b)
  colnames(x) = colnames(M)

  for(g in rownames(b)) {
    cat(g, "")
    try({
      # remove zeros implied whenever b == 0
      a = remove.zeros(M, b[g,], b.var[g,])

      # estimate just the non-zero cells
      r = lsei(A = a$M / sqrt(a$b.var), B = a$b / sqrt(a$b.var),
        G = Diagonal(sum(a$nz)), H = rep(0, sum(a$nz)), type=2)

      # only write the cells in non-zero fractions
      x[ g, a$nz ] = r$X
    })
  }

  x
}

x.pseudoinverse = { 
  x = r.corrected$r.mean %*% pseudoinverse(t(M))

  # scaling to get "average read depth / cell"
#  x = 1341 * t( t(x) / as.vector(M["all",]) )

#  x[,"P0"] = 0
  x[ is.na(x) ] = 0
  x[ x < 0 ] = 0
  x
}

test1 = function() {
  x11()
  par(mfrow=c(5,1))
#  avg.expr = apply(r1, 1, mean)
#  set.seed(0)
#  genes = sample(names(avg.expr)[avg.expr > 100], 5)
  genes = c("pha-4", "ceh-26", "pal-1", "rgs-3")
  r = unmix.lsei(M, r.corrected$r.mean[genes,], r.corrected$r.var[genes,])
  for(g in genes) {
    plot(r[g,], type="h", main=g)
  }
  r
}

# FIXME: move this elsewhere?
embryo.timeseries = read.table("git/unmix/seq/timing/embryo.timeseries.tsv.gz",
  sep="\t", header=TRUE, row.names=1, as.is=TRUE)
time.points = c(0,30,60,90,120,140,180,240,270,300,330,360,390,
  420,450,480,540,570,600,630,660,690,720)
# colnames(embryo.timeseries) = time.points
embryo.timeseries = embryo.timeseries[ rownames(r) , ]

# normalize to "ppm"
embryo.timeseries =
  t( t(embryo.timeseries) / (apply(embryo.timeseries,2,sum) / 1e6) )
g1 = intersect(rownames(embryo.timeseries), rownames(r.corrected$r.mean))

pseudoinverse.with.time = function() {
  load("git/unmix/image/time_sort_matrix.Rdata")

# scale rows of this to add up to 1
  M.t = time.sort.matrix / apply(time.sort.matrix, 1, sum)
  M.t = rbind(M, M.t)

#  x[,"P0"] = 0
  x[ is.na(x) ] = 0
  x[ x < 0 ] = 0
  x
}

# foo = test1()

