# Does unmixing, leaving out reporter genes.

library(clinfun)

source("git/unmix/unmix.r")
source("git/unmix/unmix_lsei.r")

load("git/unmix/image/cell_weights.Rdata")

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
#   unmix.function - a function which unmixes
#   b, b.var - mean and variance of expression in those fractions
# Returns: list containing
#   x - the unmixed expression predictions
#   num.reporters - number of reporters used for each
unmix.xval = function(M, unmix.function, b, b.var) {
  x = NULL
  num.reporters = c()

  # the names of genes (if any) to which fractions correspond
  fraction.gene.names = sub("_.*$", "", rownames(M))

  for(g in rownames(b)) {
cat(g, "\n")
    # remove fractions presumed to correspond to this gene
    fractions = fraction.gene.names != g

    x.p = unmix.function(M[ fractions, ], b[ g, fractions , drop=FALSE],
      b.var[ g, fractions , drop=FALSE ])
    x = rbind(x, x.p)
    num.reporters = c(num.reporters, sum(fractions))
  }

  rownames(x) = rownames(b)
  colnames(x) = colnames(M)
  names(num.reporters) = rownames(b)

  list(x = x, num.reporters = num.reporters)
}

# Runs EP for one gene.
unmix.ep.1 = function(M, b, b.var) {
  r = approx.region(M, as.vector(b), as.vector(b.var), prior.var=1e3 * max(b.var))
  a = mvnorm.2(r$m, r$v, M, as.vector(b), as.vector(b.var))

  ep.V[[ length(ep.V)+1 ]] <<- a$V  # XXX

  r$m
}

reporters = unique(sub("_.*$", "", rownames(M)))
reporters = setdiff(reporters, "all")

reporters = intersect(reporters, rownames(expr.cell))
reporters = intersect(reporters, rownames(r.corrected$m))

r.pseudo = unmix.xval(M, unmix.pseudoinverse,
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
get.processed.data = function(time, negatives, volume, purity) {
  M = NULL
  b = NULL
  b.var = NULL

  # include time?
  if (time)
    M = rbind(sort.matrix.unweighted, time.sort.matrix.unweighted)
  else
    M = sort.matrix.unweighted

  # include negatives?
  if (!negatives) {
    M = M[ grep("minus", rownames(M), invert=TRUE) , ]
  }

  # include weighting by cell volume?
  if (volume) {
    M = t( t(M) * cell.weights[ colnames(M) , "w" ] )
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
  r1 = intersect(rownames(M), colnames(b))
  M = M[r1, cells.to.include]
  M = t( t(M) / apply(M, 1, sum) )   # normalize rows to sum to 1

  list(M = M[r1, cells.to.include], b = b[,r1], b.var = b.var[,r1])
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

# Computes AUC accuracy.
auc.accuracy = function(on.off, x.prediction) {
  r = rep(NA, nrow(on.off))

  for(i in 1:nrow(on.off))
    r[i] = roc.area.test(x.prediction[i,], on.off[i,])$area
  mean(r, na.rm=TRUE)
}

# Computes accuracy by several methods.
# Args:
#   expr, M - the actual expression, and sort matrix
#   x.prediction - the expression prediction
# Returns: vector with Pearson and Spearman correlations, and AUC
compute.accuracy = function(expr, M, x.prediction) {
  g = intersect(rownames(expr), rownames(x.prediction))
  cells = intersect(colnames(expr), colnames(x.prediction))
  list( pearson = mean(diag(cor(t(expr[g,cells]), t(x.prediction[g,cells]),
      method = "pearson", use="pairwise.complete.obs"))),
    spearman = mean(diag(cor(t(expr[g,cells]), t(x.prediction[g,cells]),
      method = "spearman", use="pairwise.complete.obs"))),
    auc = auc.accuracy(M, x.prediction) )
}


# Runs several unmixing methods on the reporter genes (cross-validated),
# and computes statistics.
unmix.crossval.stats = function(method.name, unmix.function) {
  r = NULL

  for(time in c(FALSE, TRUE))
    for(negatives in c(FALSE, TRUE))
      for(volume in c(FALSE, TRUE))
        for(purity in c(FALSE, TRUE)) {
          p = get.processed.data(time, negatives, volume, purity)
          x.p = unmix.xval(p$M, unmix.function, p$b, p$b.var)
          accuracy = compute.accuracy(expr.cell[reporters, cells.to.include],
            sort.matrix.unweighted[reporters,cells.to.include] >= 0.5,
            x.p$x[reporters, cells.to.include])
          r = rbind(r, c(time=time, negatives=negatives,
            volume=volume, purity=purity,
            pearson.corr = accuracy$pearson,
            spearman.corr = accuracy$spearman,
            area.under.curve = accuracy$auc) )
        }

  cbind(method=method.name, data.frame(r))
}

crossval.accuracy.summary =
  rbind(unmix.crossval.stats("pseudoinverse", unmix.pseudoinverse),
    unmix.crossval.stats("bounded pseudoinverse", unmix.lsei),
    unmix.crossval.stats("EP", unmix.ep.1))

# unmix.crossval.stats("pseudoinverse", unmix.pseudoinverse)

write.table(crossval.accuracy.summary,
  file="git/unmix/eval/crossval_accuracy_summary.tsv",
  sep="\t", col.names=NA)


