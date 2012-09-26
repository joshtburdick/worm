# Tests effects of gating errors are on accuracy.

library("corpcor")

# source("import.r")
source("deconvolve.r")

reporters = read.csv("out/reporters/P0.csv", stringsAsFactors=FALSE)

# reporter matrix, using a fixed cutoff
m = ifelse(expr.cell >= 2000, 1, 0)
m = rbind(1, m)
rownames(m)[1] = "all"
m1 = m[ c("all", reporters$embryo.movie.name), ]

# Measures accuracy with potential gating errors included.
# Args:
#   expr - the expression matrix
#   m - the matrix indicating which cell is in which fraction
#   b.noise - parameter for the beta indicating noise
#   b.estimate - what to use as the estimate for the Beta of noise
# Returns: list with elements
#   sim.frac - simulated amount of each gene in each fraction
#   sim.expr.d - result of deconvolving, using "b.estimate" as
#     the matrix for deconvolution
sim.with.gating.errors = function(expr, m, b.noise, b.estimate) {

  # compute perturbed matrix, and find fractions
  m1 = m
  m1[ m1==0 ] = rbeta(sum(m1==0), 1, b.noise)
  sim.frac = expr %*% t(m1)

  m.d = m
  m.d[ m.d==0 ] = 1 / (1 + b.estimate)
  sim.expr.d = sim.frac %*% pseudoinverse(t(m.d))

  list(sim.frac = t(sim.frac), sim.expr.d = sim.expr.d)
}

# Computes accuracy with a variety of levels of gating errors.
# Args:
#   b.s - the parameters for the beta
# Returns: an array with dimensions:
#   b.sim - the beta used in simulating noise
#   b.assumed - the beta assumed in deconvolving
#   accuracy.measure - either "cor" or "auc"
compute.accuracy.with.gating.errors = function(b.s) {
  r = array(dim=c(length(b.s), length(b.s), 2),
        dimnames=list(b.sim=b.s, b.assumed=b.s,
          accuracy.measure = c("cor", "auc")))

  for(i in b.s)
    for(j in b.s) {
      sim = sim.with.gating.errors(expr.cell, m1, i, j)
      r[as.character(i),as.character(j),"cor"] =
        cor1(expr.cell, sim$sim.expr.d)
      r[as.character(i),as.character(j),"auc"] =
        mean(auc(ifelse(expr.cell>=2000, 1, 0), sim$sim.expr.d), na.rm=TRUE)
  }

  r
}

accuracy.with.gating.errors =
  compute.accuracy.with.gating.errors(c(3,5,7,10,100))


