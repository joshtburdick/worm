# Attempt at predicting just the mean using
# maximum likelihood.

source("git/utils.r")
source("git/unmix/ml/unmix_sort_fraction_ml_1.r")

# the sort matrix
source("git/sort_paper/unmix/sortMatrix.r")
# ??? tweak the ceh-6 & hlh-16 so that the numbers
# are 0 or 1 ? (currently one cell is at 0.5)
m1 = {
  m.plus = m.unnormalized[1:14,]
  rownames(m.plus) = paste(rownames(m.plus), "(+)")
  # XXX omitting (-) sort fractions for hlh-16 & irx-1
  m.minus = 1 - m.unnormalized[c(1:5,8:14),]
  rownames(m.minus) = paste(rownames(m.minus), "(-)")
  m1 = rbind(m.plus, m.minus, m.unnormalized[15:18,])

  # tack on "residual" 
#  resid = diag(nrow(m1))
#  colnames(resid) = paste("r", rownames(m1))
#  cbind(m1, resid)
  m1
}

# normalize m1
m1 = m1 / apply(m1, 1, sum)

# the read data
rpm = read.tsv("git/cluster/readsPerMillion.tsv")
# XXX average cnd-1 and pha-4
rpm$"cnd-1 (+)" = apply(
  rpm[,paste(c("cnd-1 12/14", "cnd-1 1/4", "cnd-1 8/19"), "(-)")],
  1, mean)
rpm$"cnd-1 (-)" = apply(
  rpm[,paste(c("cnd-1 12/14", "cnd-1 1/4", "cnd-1 8/19"), "(-)")],
  1, mean)
rpm$"pha-4 (+)" = apply(
  rpm[,paste(c("pha-4 12/9", "pha-4 5/9", "pha-4 9/1"), "(+)")],
  1, mean)
rpm$"pha-4 (-)" = apply(
  rpm[,paste(c("pha-4 12/9", "pha-4 5/9", "pha-4 9/1"), "(-)")],
  1, mean)
rpm = as.matrix(rpm[ , rownames(m1) ])

# normalize columns
rpm = t( t(rpm) / apply(rpm, 2, sum) )

# just using the genes with highest max expression
# in some sort fraction
max.expr = sort(apply(rpm, 1, max), decreasing=TRUE)
rpm = rpm[ names(max.expr)[ 1:5000 ] , ]

# tack on sum of other genes (in rpm)
x.other = 1 - apply(rpm, 2, sum)
rpm = rbind(rpm, "_other" = x.other)

# Definition of "enrichment."
log2.enrich = function(x.plus, x.minus)
  log2( 3 + x.plus ) - log2( 3 + x.minus )

# Unmixes a data set.
# Args:
#   m - the sort matrix
#   r - the reads-per-million
#   max.iters - maximum number of iterations for optimization
# Returns: the unmixed data.
unmix.ml = function(m, rpm, max.iters=200) {

  # start with "truncated pseudoinverse"
  x = pseudoinverse(m) %*% t(rpm)
  x[ x < 0 ] = 0
  rownames(x) = colnames(m)
  colnames(x) = rownames(rpm)

  r = pos.linear.solve(m, rpm, x,
    max.iters=10, eps=1e-10)

  r$X
}

# Simulates a sorting experiment.
# Args:
#   x - the predicted expression data
#   m - the sort matrix
#   f.plus, f.minus - names of positive and negative
#     fractions used for sorting
# Returns: vector of predicted read ratios, computed
#   as for the actual data.
sim.sort.ratio = function(x, m, f.plus, f.minus) {
  x.plus = as.vector( x %*% m[ f.plus , ] )
  x.minus = as.vector( x %*% m[ f.minus , ] )

  r = log2( 3 + x.plus ) - log2( 3 + x.minus )
  names(r) = rownames(x)
  r
}

# Measures accuracy of unmixing by cross-validation.
# Args:
#   rpm - the reads to use for unmixing
#   m - the sort matrix
#   f - the sort fraction to unmix
# Returns: list containing
#   x - the unmixed data
#   sim.r - the simulated enrichment
#   actual.r - the actual enrichment
unmix.ml.crossval = function(rpm, m, f) {
  # FIXME deal with cases w/o negative controls
  f.plus = paste0(f, " (+)")
  f.minus = paste0(f, " (-)")

  # remove fraction being predicted from the input data
  # (as well as doubly-sorted fractions)
  s = setdiff(rownames(m),
    c(double.sorted.fractions, f.plus, f.minus))
print(s)
#  x = t(unmix.ml(m[ s , ], rpm[ , s ]))

  x = 1e6 * t( unmix.expr.and.sort.matrix.1(m[ s , ], t(rpm[ , s ]))$x )

  # name of the fraction to use for sim. sorting
  sim.f = f
  if (f %in% single.fractions) {
    sim.f = paste(f, "(+)")
  }

  list(x = x,
    sim.r = sim.sort.ratio(x, m, f.plus, f.minus),
    actual.r = log2.enrich(rpm[ , f.plus ], rpm[ , f.minus ]))
}

unmix.ml.all = function() {
  r = NULL

  for(f in setdiff(rownames(m)[1:14], c("hlh-16", "irx-1"))) {
    write.status(f)
    a = unmix.ml.crossval(rpm, m1, f)
    r[[f]] = cor(a$sim.r, a$actual.r)
  }

  r
}

# x1 = unmix.expr.and.sort.matrix.1(m1, t(rpm))

# crossval.r = unmix.ml.crossval(rpm, m1, "mls-2")

unmix.corr = unmix.ml.all()


