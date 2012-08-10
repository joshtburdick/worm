# Does unmixing, leaving out reporter genes.

source("git/unmix/unmix.r")
source("git/unmix/unmix_lsei.r")

source("git/plot/plot_expr.r")

# XXX for storing full covariance matrices. This is a hack.
ep.V = NULL

# tack on newer mls-2 expression data
expr.cell = {
  a = expr.cell
  load("git/image/expr.cell.Rdata")
  a = rbind(a, expr.cell)
}

# Computes accuracy of an unmixing method on a set of genes.
# Args:
#   M - a cell-sorting matrix (assumes that if a fraction corresponds
#     to a gene, then chopping off everything after "_" will give that
#     gene name. see code if that's confusing.)
#   expr - actual expression of a set of genes
#   unmix.function - a function which unmixes
#   b, b.var - mean and variance of expression in those fractions
# Returns: list containing
#   x - the unmixed expression predictions
#   num.reporters - number of reporters used for each
unmix.xval = function(M, expr, unmix.function, b, b.var) {
  x = NULL
  num.reporters = c()

  # the names of genes (if any) to which fractions correspond
  fraction.gene.names = sub("_.*$", "", rownames(M))

  for(g in rownames(expr)) {
cat(g, "\n")
    # remove fractions presumed to correspond to this gene
    fractions = fraction.gene.names != g

    x.p = unmix.function(M[ fractions, ], b[ g, fractions , drop=FALSE],
      b.var[ g, fractions , drop=FALSE ])
    x = rbind(x, x.p)
    num.reporters = c(num.reporters, sum(fractions))

#    r = rbind(r, c(cor = cor(expr[g,], x[1,], use="complete.obs"),
#      num.reporters = sum(fractions) ))
  }

  rownames(x) = rownames(expr)
  colnames(x) = colnames(M)
  names(num.reporters) = rownames(expr)

  list(x = x, num.reporters = num.reporters)
}

# Runs EP for one gene.
unmix.ep.1 = function(M, b, b.var) {
  r = approx.region(M, as.vector(b), as.vector(b.var), prior.var=1e6 * max(b.var))
  a = mvnorm.2(r$m, r$v, M, as.vector(b), as.vector(b.var))

  ep.V[[ length(ep.V)+1 ]] <<- a$V  # XXX

  r$m
}

reporters = unique(sub("_.*$", "", rownames(M)))
reporters = setdiff(reporters, "all")

reporters = intersect(reporters, rownames(expr.cell))
reporters = intersect(reporters, rownames(r.corrected$m))

r.pseudo = unmix.xval(M, expr.cell[reporters,colnames(M)], unmix.pseudoinverse,
  r.corrected$m[reporters,], r.corrected$v[reporters,])

# Gets a particular version of the cell-sorting matrix, and corresponding
# measurements in fractions, with and without various corrections.
# Args:
#   time - whether or not to include time
#   negatives - whether or not to include the negative fractions
#   volume - whether or not to include cell-volume correction
#   purity - whether or not to include the sort purity correction
# Returns: list containing:
#   m - a cell-sorting matrix
#   b, b.var - the corresponding total expression (mean and variation)
get.input.data = function(time, negatives, volume, purity) {
  m = NULL
  b = NULL
  b.var = NULL

  # include time?
  if (time)
    m = rbind(sort.matrix, time.sort.matrix)
  else
    m = sort.matrix

  # include negatives?
  if (!negatives) {
    m = m[ grep("minus", rownames(m), invert=TRUE) , ]
  }

  # include weighting by cell volume?
  if (volume) {
    m = t( t(m) * cell.weights[ colnames(m) , "w" ] )
  }

  # include sort purity correction?
  if (purity) {
    b = cbind(r.corrected$m[reporters,], r.ts$m[reporters,])
    b.var = cbind(r.corrected$v[reporters,], r.ts$v[reporters,])
  }
  else {
    b = cbind(r[reporters,], r.ts$m[reporters,])
    b.var = b   # assuming the variance is the same as the mean
  }

  stopifnot(all(colnames(b) == colnames(b.var)))
  r1 = intersect(rownames(m), colnames(b))

  list(m = m[r1, cells.to.include], b = b[,r1], b.var = b.var[,r1])
}

if (FALSE) {
  r.lsei = unmix.xval(M, expr.cell[reporters,colnames(M)], unmix.lsei,
  r.corrected$m[reporters,], r.corrected$v[reporters,])
}

if (FALSE) {
  ep.V <<- list()
  r.ep = unmix.xval(M, expr.cell[reporters,colnames(M)], unmix.ep.1,
    r.corrected$m[reporters,], r.corrected$v[reporters,])
  names(ep.V) = reporters
}

# Computes accuracy by several methods.
# Args:
#   expr, m - the actual expression, and sort matrix
#   x.prediction - the expression prediction
# Returns: vector with Pearson and Spearman correlations, and AUC
compute.accuracy = function(expr, m, x.prediction) {
  g = intersect(rownames(x1), rownames(x2))
  cells = intersect(colnames(x1), colnames(x2))
  c( pearson = mean(diag(cor(t(x1[g,cells]), t(x2[g,cells]),
      method = pearson, use="pairwise.complete.obs"))),
    spearman = mean(diag(cor(t(x1[g,cells]), t(x2[g,cells]),
      method = separman, use="pairwise.complete.obs"))),
    auc = 0 )     # FIXME
}

