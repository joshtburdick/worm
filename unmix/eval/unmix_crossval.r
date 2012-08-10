# Does unmixing, leaving out reporter genes.

source("git/unmix/unmix.r")
source("git/unmix/unmix_lsei.r")

source("git/plot/plot_expr.r")

# XXX for storing full covariance matrices
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
# Finds correlation with the reporters.
compute.correlation = function(x1, x2) {
  g = intersect(rownames(x1), rownames(x2))
  cells = intersect(colnames(x1), colnames(x2))
  mean(diag(cor(t(x1[g,cells]), t(x2[g,cells]), use="pairwise.complete.obs")))
}


