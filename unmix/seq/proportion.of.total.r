# Converts from read depth to "proportion of expression in a fraction."

# mean and variance of a beta distribution
beta.mean = function(a, b) a / (a + b)
beta.variance = function(a, b) (a*b) / ( (a+b)^2 * (a+b+1) )

# Estimates proportion of expression in given fractions,
# based on read data.
# Args:
#   r - the read coverage in each fraction. The fraction
#     containing (presumably) all cells should be called "all",
#     and negative fractions should have a suffix of "_minus".
#   fraction.of.cells.sorted - proportion of cells included
#     in each fraction (according to imaging data)
# Returns: a list of two matrices "m" and "v", giving the
#   mean and variance of the total expression in a fraction
#   of each gene (thus, rescaled to the range [0,1].)
read.depth.to.fraction = function(r, fraction.of.cells.sorted) {
  prop.expr = r

  # For fractions where we don't have negatives, estimate
  # "proportion of total expression" relative to control.
  get.proportion.pos.only = function(prop.expr) {

    # proportion of expression in the positive fraction
    a = t( fraction.of.cells.sorted[-1] * t(prop.expr[,-1] / prop.expr[,1]) )
    a[ is.na(a) ] = 0
    a[ a >= 1 ] = 1

    # total reads relevant to each gene
    tr = (r1[,-1] + r1[,1])

    list(m = beta.mean(a*tr + 1, (1-a)*tr + 1),
      v = beta.variance(a*tr + 1, (1-a)*tr + 1))
  }

  # For fractions where we do have negatives, use the
  # ratio to estimate a proportion.
  get.proportion.pos.neg = function(prop.expr, gene) {
    pos = gene
    neg = paste(gene, "_minus", sep="")

    # proportion of reads in the positive fraction
    p = fraction.of.cells.sorted[pos]
    a = p * prop.expr[,pos] / (p * prop.expr[,pos] + (1-p) * prop.expr[,neg] )
    a[ is.na(a) ] = 0

    # total reads for this gene
    tr = (r1[,pos] + r1[,neg])

    list(m = beta.mean(a*tr + 1, (1-a)*tr + 1),
      v = beta.variance(a*tr + 1, (1-a)*tr + 1))
  }

  # find cases in which we have negative fractions
  pos.neg.fractions = sub("_minus", "", grep("minus", colnames(prop.expr)))
  pos.only.fractions = setdiff(colnames(r), pos.neg.fractions)
  pos.only.fractions = setdiff(pos.only.fractions, "all")

  a = get.proportion.pos.only(r[,c("all", pos.only.fractions)])
  for(g in pos.neg.fractions) {
    a1 = get.proportion.pos.neg(prop.expr, g)
    a$m <- cbind(a$m, a1$m)
    a$v <- cbind(a$v, a1$v)
  }

  a
}

