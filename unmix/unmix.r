# Gets all the data sets for unmixing.

library(corpcor)
library(limSolve)

load("git/unmix/image/sort_matrix.Rdata")

# the expression data as read depth, with correction for sorting purity
source("git/unmix/seq/sortPurityCorrection.r")
names(r.corrected) = c("m", "v")

# limit sort matrix to only include things up to a certain time
t.cutoff = 420
cell.time.on.off =
  read.csv(gzfile("data/image/cellTimeOnOff.csv.gz"),
    row.names=1, as.is=TRUE)

# which cells to include
# XXX note that this is not in lineage order
cells.to.include =
  rownames(cell.time.on.off)[ cell.time.on.off$on < t.cutoff ]
M = sort.matrix[ colnames(r.corrected$m), cells.to.include ]

# scale rows of this to add up to 1
M = M / apply(M, 1, sum)

# limit to cases in which we have measurements
M = M[ colnames(r.corrected$m) , ]

# sorting matrix including temporal measurements
load("git/unmix/image/time_sort_matrix.Rdata")
# only include times <= 420, and that set of cells
M.t = time.sort.matrix[ c(1:13) , cells.to.include ]
M.t = M.t / apply(M.t, 1, sum)

# the sort matrix including the FACS data and the timeseries
M.facs.and.ts = rbind(M, M.t)

# the timeseries of embryonic expression data
embryo.timeseries = as.matrix(
  read.table("git/unmix/seq/timing/embryo.timeseries.tsv.gz",
    sep="\t", header=TRUE, row.names=1, as.is=TRUE))
time.points = paste("t.", c(0,30,60,90,120,140,180,240,270,300,330,360,390,
  420,450,480,540,570,600,630,660,690,720), sep="")
colnames(embryo.timeseries) = time.points

# normalize to "ppm"
embryo.timeseries =
  t( t(embryo.timeseries) / (apply(embryo.timeseries,2,sum) / 1e6) )

g1 = intersect(rownames(embryo.timeseries), rownames(r.corrected$m))

# again, we assume that the variance in this is about the same as
# expression on the ppm scale
r.ts = list(m = embryo.timeseries[ g1, rownames(M.t) ],
  v = embryo.timeseries[g1 , rownames(M.t) ])
r.corrected = list(m = r.corrected$m[g1,], v = r.corrected$v[g1,])

# get list of genes to unmix
load("R/unmix/comp_paper/expr.cell.Rdata")
enriched.fraction =
  read.table("R/unmix/sort_paper/unmix/fraction/enriched.fraction.tsv",
  sep="\t", header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE)
gene.list = c(rownames(expr.cell), rownames(enriched.fraction))
gene.list = unique(intersect(gene.list, rownames(r.corrected$m)))

# Removes cells in any fractions having zero expression.
# Args:
#   M - the sort matrix
#   b - expression in each fraction
# Returns: a list with components
#   M, b, b.var - the corresponding variables, with zeros
#     removed from b (and M and b.var; M will also have
#     columns corresponding to zero fractions removed)
#   nz - which columns were nonzero
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
      # remove zeros implied whenever b == 0
#      a = remove.zeros(M, b[g,], b.var[g,])

      # compute prior, with the fuzzy constraint observed
      p = prior.var * max(b.var[g,], na.rm=TRUE)
      mv = mvnorm.2(rep(0, ncol(M)), rep(p, ncol(M)),
        M, b[g,], b.var[g,])
print(mv$V[1:5,1:5])
print(mv$m[1:10])
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

# Unmix using the pseudoinverse.
# Args:
#   M - the cell-sorting matrix
#   b - the expression in each fraction
# Returns: matrix of predictions (truncated to be positive.)
unmix.pseudoinverse = function(M, b) {
  r = b %*% pseudoinverse(t(M))
  r[ r < 0 ] = 0
  r
}

x.pseudo = unmix.pseudoinverse(M, r.corrected$m)
x.pseudo.time = unmix.pseudoinverse(M.facs.and.ts,
  cbind(r.corrected$m, r.ts$m))

unmix.lsei.1 = function() {
  x11()
  par(mfrow=c(5,1))
#  avg.expr = apply(r1, 1, mean)
#  set.seed(0)
#  genes = sample(names(avg.expr)[avg.expr > 100], 5)
  genes = c("pha-4", "ceh-26", "pal-1", "rgs-3")
  r = unmix.lsei(M, r.corrected$m[genes,], r.corrected$v[genes,], prior.var=1e2)
#  for(g in genes) {
#    plot(r[g,], type="h", main=g)
#  }
  r
}


