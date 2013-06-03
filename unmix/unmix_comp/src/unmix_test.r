# Runs an unmixing method on several datasets,
# using various numbers of reporters.

library("limSolve")

expr.cell = as.matrix(read.table("../data/exprCell.tsv",
  sep="\t", header=TRUE, row.names=1, as.is=TRUE))

load("sim.data/sim.cell.cor.Rdata")
load("sim.data/sim.expr.cell.Rdata")
load("sim.data/sim.expr.sym.Rdata")

source("expression.on.off.r")

# cell-sorting matrix
m.cell = rbind(all=1, expression.on.off(expr.cell))

reporters = read.table("reporters.tsv",
  row.names=1, header=TRUE, as.is=TRUE)

num.reporters = c(10,20,30,50,75,100)

# Runs an unmixing function with different numbers of reporters.
# When predicting a given gene, avoids using it as a reporter.
run.unmix.1 = function(x, m, unmix.f, reporter.list, nr) {
  unmix.result = list(
    x.predicted = matrix(nrow=nrow(x), ncol=ncol(x)),
    r = list()
  )

  dimnames(unmix.result$x.predicted) = dimnames(x)

  for(g in rownames(x)) {
cat(g, "\n")
      # set of reporters, with another one added in if g is present
    r1 = if (g %in% reporter.list[1:nr])
      setdiff(reporter.list[1:(nr+1)], g)
    else
      reporter.list[1:nr]
# cat("r1 =", r1, "\n")

    m1 = m[ c("all", r1), ]
    x.fraction = x[g,] %*% t(m1)

    # this is wrapped in a "try" block, in case unmixing fails
    try({
      r = unmix.f(m1, as.vector(x.fraction))
      unmix.result$x.predicted[g,] = as.vector(r$x)
      unmix.result$r[[g]] = r
    })
  }

  rownames(unmix.result$x.predicted) = rownames(x)

  unmix.result
}

# Does deconvolution on several data sets.
# Args:
#   unmix.f - a function which, given a cell-fraction matrix m, and a
#     set of fractions x.f, returns a vector of unmixed expression
#   output.dir - where to write output
run.unmix = function(unmix.f, output.dir) {
  system(paste("mkdir -p ", output.dir))

  run.1 = function(expr, name) {
    for(nr in num.reporters) {
      unmix.result = NULL
      st = system.time( unmix.result <-
        run.unmix.1(expr, m.cell, unmix.f, reporters$picked, nr) )
      unmix.result$system.time = st
      save(unmix.result,
        file=paste(output.dir, name, ".", nr, ".Rdata", sep=""))
    }
  }

  run.1(expr.cell, "expr.cell")
  run.1(sim.expr.cell, "sim.expr.cell")
  run.1(sim.cell.cor, "sim.cell.cor")
  run.1(sim.expr.sym, "sim.expr.sym")
}

