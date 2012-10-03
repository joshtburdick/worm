# Tests effects of missing cells on accuracy.

library("corpcor")

# source("import.r")
# source("deconvolve.r")

source("~/gcb/R/unmix/eval.r")

load("~/gcb/git/unmix/unmix_comp/data/tree_utils.Rdata")

wd = getwd()
setwd("~/gcb/git/unmix/unmix_comp/src/")
source("unmix_test.r")
setwd(wd)

num.cells = apply(cell.lineage.matrix, 1, sum)

# Creates a random reporter matrix.
# Args:
#   num.reporters - the number of reporters
#   p.on - probability of a cell being "on"
# Returns: a random sort matrix
random.sort.matrix = function(num.reporters, p.on) {
  a = matrix(rbinom(1341 * num.reporters, 1, p.on),
    nrow=num.reporters, ncol=1341)
  a[1,] = 1
  colnames(a) = lin.node.names
  a
}

m1 = m.cell[ rownames(reporters[1:30,]) , ]

# Computes accuracy using random reporters
# Args:
#   x - the expression data to use for testing
#   num.reporters - vector of number of reporters to test
#   p.on - vector of probability on to test
#   num.reps - number of different matrices to test
# Returns: data frame of results
compute.random.sort.accuracy =
    function(x, num.reporters, p.on, num.reps) {
  r = NULL

  for(nr in num.reporters)
    for(p in p.on)
      for(i in 1:num.reps) {
cat(nr, p, "\n")
        m.random = random.sort.matrix(nr, p)
        x.d = t( pseudoinverse(m.random) %*% m.random %*% t(x) )
        r = rbind(r, c(num.reporters = nr, p.on = p,
          cor = cor1(expr.cell, x.d),
          auc = mean(auc(expr.cell >= 2000, x.d))))
      }

  r
}

random.matrix.accuracy = compute.random.sort.accuracy(
  expr.cell, 10 * c(1:7), c(0.1, 0.5, 0.9), 5)
write.table(random.matrix.accuracy,
  file="git/unmix/random_reporter_accuracy.tsv", sep="\t",
  row.names=FALSE, col.names=TRUE)

