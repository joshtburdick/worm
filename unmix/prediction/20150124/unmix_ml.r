# Attempt at predicting just the mean using
# maximum likelihood.

source("git/utils.r")
source("git/unmix/ml/pos_linear_solve.r")

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

# restrict to genes with at least some amount of expression)
expressed = which( apply(rpm, 1, max) >= 5 )

# only include a random subset of these?
# set.seed(42)
# rpm = rpm[ sample(expressed, 3000) , ]

# tack on sum of other genes (in rpm)
# x.other = 1e6 - apply(rpm, 2, sum)
# rpm = rbind(rpm, "_other" = x.other)

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
sim.sort.ratio = function(x, m) {
  m.plus = m / sum(m)
  m.minus = (1-m) / sum(1-m)

  x.plus = x %*% as.matrix(as.vector(m.plus))
  x.minus = x %*% as.matrix(as.vector(m.minus))

  log2( 3 + x.plus ) - log2( 3 + x.minus )
}

# Measures accuracy of unmixing by cross-validation.
# Args:
#   rpm - the reads to use for unmixing
#   f - the sort fraction to unmix
# Returns: list containing
#   x - the unmixed data
#   sim.r - the simulated 
unmix.ml.crossval = function(rpm, f) {

  # remove fraction being predicted from the input data
  # (as well as doubly-sorted fractions)
  s = setdiff(rownames(m1),
    c(double.sorted.fractions, paste(f, c("(+)", "(-)"))))
# print(s)
  x = unmix.ml(m1[ s , ], rpm[ , s ])

  # name of the fraction to use for sim. sorting
  sim.f = f
  if (f %in% single.fractions) {
    sim.f = paste(f, "(+)")
  }

  list(x = x,
    sim.r = sim.sort.ratio(x[,1:1341], m1[sim.f, 1:1341]))
}



# x = unmix.ml(m1, rpm)

# crossval.r = unmix.ml.crossval(rpm, "pros-1")

