# Picking reporters to maximize accuracy.

library(corpcor)

load("../data/tree_utils.Rdata")

expr.cell = as.matrix(read.table("../data/exprCell.tsv",
  sep="\t", header=TRUE, row.names=1, as.is=TRUE))

load("sim.data/sim.cell.cor.Rdata")
load("sim.data/sim.expr.cell.Rdata")
load("sim.data/sim.expr.sym.Rdata")

source("expression.on.off.r")
# cell-sorting matrix
m.cell = rbind(all=1, expression.on.off(expr.cell))

# Utility to standardize the rows of a matrix.
standardize.rows = function(x)
  (x - apply(x,1,mean)) / apply(x,1,sd)

# Computes correlation of each row of two matrices.
# The rows of the matrices must be standardized.
row.cor = function(x, y) {
  n = dim(x)[2]
  apply(x*y, 1, sum) / (n-1)
}

# Computes AUC.
compute.auc = function(x.on.off, x.prediction) {
  x.on.off = ifelse(as.matrix(x.on.off) >= 0.5, 1, 0)

  n = dim(x.on.off)[1]
  auc = rep(NA, n)

  for(i in 1:n) {
    if (sd(x.on.off[i,]) > 0 && sum(is.na(x.prediction[i,])==0)) {
#        cat(i, "")
      a = roc.area.test(x.prediction[i,], x.on.off[i,])
      auc[i] = as.numeric(a$area)
    }
  }

  mean(auc, na.rm=TRUE)
}

# Computes the accuracy of some set of reporters.
# For now, "accuracy" is "correlation in some cells".
# Args:
#   m - the matrix of which cell is in which fraction
#   x - the expression dataset to use in measuring accuracy
#   focus - which cells (columns) to include in correlation
# Returns: a function from lists of reporters to accuracy.
reporter.accuracy = function(m, x, focus) {
  x.s = standardize.rows(x[,focus])

  function(reporters) {
    m1 = as.matrix(m[,reporters])
    x.fraction = x %*% m1
    x.u = x.fraction %*% pseudoinverse(m1)
    x.u.s = standardize.rows(x.u[,focus])
    mean(row.cor(x.s, x.u.s))
# XXX trying using AUC criterion
#    compute.auc(sim.expr.cell / 10, x.u.s)
  }
}

# Picks a set of things to maximize some objective function.
# Args:
#   f - the function to maximize
#   xs - list of things to pick from
#   num.to.pick - number of things to pick
# Returns: a data frame with two columns
pick.greedily = function(f, xs, num.to.pick) {
  picked = c()
  f.val = c()

  for(iter in 1:num.to.pick) {

    # evaluate f with each thing added, and pick maximum
#    f1 = function(x) f(c(picked, x))
#   verbose variant
    f1 = function(x) { r = f(c(picked, x)); cat(x,r," "); r }
    r = sapply(xs, f1)
    i.max = which.max(r)

    # store which was picked, and the function value including it
    picked = c(picked, xs[i.max])
    f.val = c(f.val, r[i.max])

    cat(xs[i.max], r[i.max], "\n")

    # remove the one that was picked
    xs = setdiff(xs, xs[i.max])
  }

  data.frame(picked = picked, f.val = f.val)
}

accuracy.1 = function(r)
  reporter.accuracy(t(m.cell), sim.expr.cell, lin.node.names)(c("all", r))

reporters = pick.greedily(accuracy.1, rownames(m.cell), 110)

write.table(reporters, file="reporters.tsv",
  sep="\t", row.names=TRUE, col.names=NA)

