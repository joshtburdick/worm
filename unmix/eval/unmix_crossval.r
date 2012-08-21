# Does unmixing, leaving out reporter genes.

library(clinfun)

source("git/unmix/unmix.r")
source("git/unmix/unmix_lsei.r")
source("git/unmix/seq/proportion.of.total.r")

load("git/unmix/image/cell_weights.Rdata")

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
# cat(g, "\n")
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
  a = mvnorm.2(r$t[,"m"], r$t[,"v"], M, as.vector(b), as.vector(b.var))

#  ep.V[[ length(ep.V)+1 ]] <<- a$V        # XXX

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
#   norm.to.unit - normalize to numbers adding up to 1 ?
# Returns: list containing:
#   m - a cell-sorting matrix
#   b, b.var - the corresponding total expression (mean and variation)
get.processed.data = function(time, negatives, volume, purity, norm.to.unit) {
  M = NULL
  b = NULL
  b.var = NULL

#  cat(time,negatives,volume,purity,norm.to.unit,"\n")

  # include time?
  if (time)
    M = rbind(sort.matrix.unweighted, time.sort.matrix.unweighted)
  else
    M = sort.matrix.unweighted

  # XXX this seems to work slightly better, although I have no idea why.
  if (time)
    M = rbind(sort.matrix, time.sort.matrix)
  else
    M = sort.matrix

  M = M[,cells.to.include]

  # include negatives?
  if (!negatives) {
    M = M[ grep("minus", rownames(M), invert=TRUE) , ]
  }

  # include weighting by cell volume?
  if (volume) {
    M = t( t(M) * cell.weights[ colnames(M) , "w" ] )
  }

  # include sort purity correction?
  if (purity == "pos.neg") {
    b = cbind(r.corrected$m[reporters,], r.ts$m[reporters,])
    b.var = cbind(r.corrected$v[reporters,], r.ts$v[reporters,])
  }
  else if (purity == "unsorted.singlets") {
    load("git/unmix/seq/FACS/r.ungated.corr.Rdata")
    b = cbind(r.ungated.corr[reporters,], r.ts$m[reporters,])
    b.var = b
  }
  else {
    b = cbind(r[reporters,], r.ts$m[reporters,])
    b.var = b   # assuming the variance is the same as the mean
  }

  stopifnot(all(colnames(b) == colnames(b.var)))
  r1 = intersect(rownames(M), colnames(b))
cat(r1,"\n")
  M = M[r1,cells.to.include]
  b = b[,r1]
  b.var = b.var[,r1]

  # normalize so that we're using numbers which add up to 1 ?
  if (norm.to.unit) {
    M = M / sum(M["all",])
    fraction.of.cells.sorted = apply(M, 1, sum)
    p = read.depth.to.fraction(b, fraction.of.cells.sorted)
    b = p$m
    b.var = p$v

    # if we're normalizing in this way, omit the negatives (if any)
    r1 = grep("minus", r1, invert=TRUE, value=TRUE)
  }
  else {
#    M = t( t(M) / apply(M, 1, sum) )   # normalize rows to sum to 1
    M = M / apply(M, 1, sum)
  }

  list(M = M[r1, cells.to.include], b = b[,r1], b.var = b.var[,r1])
#  list(M = M, b = b, b.var = b.var)
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

# Runs an unmixing method on the reporter genes (with cross-validation),
# including several different options for pre-processing the data,
# and saves the result in a file.
unmix.crossval = function(unmix.function, output.file) {
  expr.prediction = list()

  for(time in c(FALSE, TRUE))
    for(negatives in c(FALSE, TRUE))
      for(volume in c(FALSE, TRUE))
        for(purity in c("pos.neg", "unsorted.singlets", " "))
          for(norm.to.unit in c(FALSE, TRUE)) {
            p = get.processed.data(time, negatives, volume, purity,
              norm.to.unit)
            x = unmix.xval(p$M, unmix.function, p$b, p$b.var)
            r = list(time=time, negatives=negatives, volume=volume,
              purity=purity, norm.to.unit=norm.to.unit,
              x.p.mean = x$x)
            expr.prediction <- c(expr.prediction, list(r))

            # FIXME: save variance, if it's present
          }
  save(expr.prediction, file=output.file)
  expr.prediction
}

# Computes accuracy for a prediction.
# Args:
#   predictions - a list of predictions (as from unmix.crossval.stats)
#   expr - the actual expression (e.g., from imaging)
#   M - cell-sorting matrix: 0-1 matrix indicating if a gene
#     is "on" or "off" in a given cell (for AUC)
# Returns: vector with Pearson, Spearman, and AUC
compute.accuracy = function(expr.prediction, expr, M) {
  accuracy = compute.accuracy(expr.cell[reporters, cells.to.include],
    sort.matrix.unweighted[reporters,cells.to.include] >= 0.5,
    x.p$x[reporters, cells.to.include])
}

if (FALSE) {

  foo = unmix.crossval(unmix.pseudoinverse,
    "git/unmix/eval/pseudoinverse_xval.Rdata")
  unmix.crossval(unmix.ep.1,
    "git/unmix/eval/ep_xval.Rdata")
#  unmix.crossval(unmix.lsei,
#    "git/unmix/eval/bound_pseudoinverse_xval.Rdata")

# write.table(crossval.accuracy.summary,
#   file="git/unmix/eval/crossval_accuracy_summary.tsv",
#   sep="\t", col.names=NA)
}


