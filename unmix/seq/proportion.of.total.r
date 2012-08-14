# Converts from read depth to "proportion of expression in a fraction."

# mean and variance of a beta distribution
beta.mean = function(a, b) a / (a + b)
beta.variance = function(a, b) (a*b) / ( (a+b)^2 * (a+b+1) )

# Estimates proportion of expression in given fractions,
# based on read data.
# Args:
#   r - the read coverage in each fraction. The fraction
#     containing (presumably) all cells should be called "all",
#     (and this should be the first column).
#     Negative fractions should have a suffix of "_minus".
#   fraction.of.cells.sorted - proportion of cells included
#     in each fraction (according to imaging data)
# Returns: a list of two matrices "m" and "v", giving the
#   mean and variance of the total expression in a fraction
#   of each gene (thus, rescaled to the range [0,1].)
read.depth.to.fraction = function(r, fraction.of.cells.sorted) {

  # For fractions where we don't have negatives, estimate
  # "proportion of total expression" relative to control.
  get.proportion.pos.only = function(r) {
    r1 = r / r[,"all"]
    fraction.of.cells.sorted = fraction.of.cells.sorted[ colnames(r1) ]

    # proportion of expression in the positive fraction
    a = t( fraction.of.cells.sorted[-1] * t(r1[,-1] / r1[,1]) )
    a[ is.na(a) ] = 0
    a[ a >= 1 ] = 1

    # total reads relevant to each gene
    tr = (r[,-1] + r[,1])

    list(m = beta.mean(a*tr + 1, (1-a)*tr + 1),
      v = beta.variance(a*tr + 1, (1-a)*tr + 1))
  }

  # For fractions where we do have negatives, use the
  # ratio to estimate a proportion.
  get.proportion.pos.neg = function(r, gene) {
    r1 = r / r[,"all"]

    pos = gene
    neg = paste(gene, "_minus", sep="")

    # proportion of reads in the positive fraction
    p = fraction.of.cells.sorted[pos]
    a = p * r1[,pos] / (p * r1[,pos] + (1-p) * r1[,neg] )
    a[ is.na(a) ] = 0

    # total reads for this gene
    tr = (r[,pos] + r[,neg])

    list(m = beta.mean(a*tr + 1, (1-a)*tr + 1),
      v = beta.variance(a*tr + 1, (1-a)*tr + 1))
  }

  # find cases in which we have negative fractions
  neg.fractions =
    sub("_minus", "", grep("minus", colnames(r), value=TRUE))
  pos.fractions = grep("minus", colnames(r),
    value=TRUE, invert=TRUE)
  pos.fractions = setdiff(pos.fractions, neg.fractions)
  pos.fractions = setdiff(pos.fractions, "all")

  # first, get positive fractions
  a = get.proportion.pos.only(r[,c("all", pos.fractions)])

  # tack on cases with negative fractions
  for(gene in neg.fractions) {
    a1 = get.proportion.pos.neg(r, gene)
    a$m <- cbind(a$m, a1$m)
    a$v <- cbind(a$v, a1$v)
    colnames(a$m)[ncol(a$m)] = gene
    colnames(a$v)[ncol(a$v)] = gene
  }

  a
}

