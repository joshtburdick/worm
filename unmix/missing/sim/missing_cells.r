# Tests effects of missing cells on accuracy.

library("corpcor")

# source("import.r")
# source("deconvolve.r")

source("~/gcb/R/unmix/eval.r")

load("~/gcb/git/unmix/unmix_comp/data/tree_utils.Rdata")

wd = getwd()
setwd("~/gcb/git/unmix/unmix_comp/src/")
source("unmix_test.r")
source("pseudoinverse/unmix.r")
setwd(wd)

num.cells = apply(cell.lineage.matrix, 1, sum)

m1 = m.cell[ c("all", rownames(reporters[1:31,])) , ]

# Runs an unmixing function with different numbers of reporters.
# When predicting a given gene, avoids using it as a reporter.
# This version also simulates data using a perturbed matrix.
run.unmix.perturbed.matrix = function(x, m, m.perturbed, unmix.f, reporter.list, nr) {
  x.predicted = matrix(nrow=nrow(x), ncol=ncol(x))
  dimnames(x.predicted) = dimnames(x)

  for(g in rownames(x)) {
 cat(g, "")

    # set of reporters, with another one added in if g is present
    r1 = if (g %in% reporter.list[1:nr])
      setdiff(reporter.list[1:(nr+1)], g)
    else
      reporter.list[1:nr]
# cat("r1 =", r1, "\n")

    m1 = m[ c("all", r1), ]
    m.perturbed.1 = m.perturbed[ c("all", r1), ]
    x.fraction = x[g,] %*% t(m.perturbed.1)

    r = unmix.f(m1, as.vector(x.fraction))
    x.predicted[g,] = as.vector(r$x)
  }

  rownames(x.predicted) = rownames(x)

  x.predicted
}

# Picks some missing cells, by picking some random lineages.
# Args: see sim.with.missing.cells()
# Returns: a Boolean vector
pick.missing.cells = function(num.lineages.missing, min.cells, max.cells) {
  lin1 = names(num.cells[ num.cells >= min.cells & num.cells <= max.cells ])
  missing.lineages = sample(lin1, num.lineages.missing)
  missing.cells = apply(cell.lineage.matrix[ missing.lineages , , drop=FALSE ], 2, sum) > 0

  missing.cells
}

# Constructs a cell-sorting matrix, with some cells missing.
# Args:
#   expr - the expression data
#   m - the cell-sorting matrix, without lineages missing
#   num.lineages.missing - number of lineages to omit
#   min.cells, max.cells - minimum and maximum number of cells to omit
#   missing.in.all.fractions - if TRUE (default), this removes cells
#     from all sort markers. if FALSE, removes cells from different fractions separately.
# Returns: list with elements:
#   m - the perturbed cell-sorting matrix
#   missing - vector of which cells were missing
sim.with.missing.cells = function(expr, m,
    num.lineages.missing, min.cells, max.cells, missing.in.all.fractions=TRUE) {
  m1 = m
  num.missing = 0
  if (missing.in.all.fractions) {
    # pick some cells to be missing, and zero them out in all fractions
    missing.cells = pick.missing.cells(num.lineages.missing, min.cells, max.cells)
    m1[ , missing.cells ] = 0
    num.missing = sum(missing.cells)
  }
  else {
    # pick some cells in individual fractions to be missing
    num.missing = 0
    mask = 1 + 0 * m1
    for(i in 1:num.lineages.missing) {
      missing.cells = pick.missing.cells(1, min.cells, max.cells)
      num.missing = num.missing + sum(missing.cells)
      mask[ sample(1:nrow(m1), 1), missing.cells ] = 0
    }
    # m1 = m1 * mask     this only zeros out some entries
    m1 = m1 * mask + (1 - m1) * (1 - mask)   # this "flips" the selected entries
    num.missing = sum( mask==0 )
  }

  # compute fractions (with the original matrix, and the perturbed matrix)
  sim.frac = expr %*% t(m)
  sim.frac.perturbed = expr %*% t(m1)

  # do deconvolving (with correct and incorrect sort matrix)  
  x.d.correct.m = run.unmix.perturbed.matrix(expr, m, m, unmix.pseudoinverse,
    rownames(m), 30)

  # XXX this is poorly named; the matrix is only "incorrect" because the
  # underlying fractions have been perturbed.
  x.d.incorrect.m = run.unmix.perturbed.matrix(expr, m, m1, unmix.pseudoinverse,
    rownames(m), 30)

  list(x.d.incorrect.m = x.d.incorrect.m, x.d.correct.m = x.d.correct.m,
    num.missing = num.missing,
    cor.incorrect.m = cor1(expr.cell, x.d.incorrect.m),
    auc.incorrect.m = mean(auc(ifelse(expr.cell>=2000, 1, 0), x.d.incorrect.m), na.rm=TRUE),
    cor.correct.m = cor1(expr.cell, x.d.correct.m),
    auc.correct.m = mean(auc(ifelse(expr.cell>=2000, 1, 0), x.d.correct.m), na.rm=TRUE))
}

# Simulates with a variety of different numbers of missing cells.
# Returns: data.frame with columns:
#   num.lineages.missing, min.cells, max.cells - used to pick missing cells
#   num.missing - number of cells which were missing
#   cor.incorrect.m, auc.incorrect.m - accuracy when using the incorrect sort matrix
#   cor.correct.m, auc.correct.m - accuracy when using the correct sort matrix
sim.with.various.params = function() {
  r = NULL

#  for(num.lineages.missing in c(1:9))
#    for(min.cells in c(1, 20, 40)) {
  for(num.lineages.missing in c(1:20))
    for(min.cells in c(1, 15, 30)) {
cat(num.lineages.missing, min.cells, "\n")
      max.cells = min.cells + 10
      a = sim.with.missing.cells(expr.cell, m1,
        num.lineages.missing, min.cells, max.cells, TRUE)
      r = rbind(r, c(num.lineages.missing = num.lineages.missing,
        min.cells = min.cells, max.cells = max.cells, num.missing = a$num.missing,
        cor.incorrect.m = a$cor.incorrect.m, auc.incorrect.m = a$auc.incorrect.m,
        cor.correct.m = a$cor.correct.m, auc.correct.m = a$auc.correct.m))
  }

  r
}

# As above, but removes different cells from different fractions.
sim.removing.from.individual.fractions = function() {
  r = NULL

  for(num.lineages.missing in c(1:70))
    for(min.cells in c(20, 40, 60)) {
cat(num.lineages.missing, min.cells, "\n")
      max.cells = min.cells + 20
      a = sim.with.missing.cells(expr.cell, m1,
        num.lineages.missing, min.cells, max.cells, FALSE)
      r = rbind(r, c(num.lineages.missing = num.lineages.missing,
        min.cells = min.cells, max.cells = max.cells, num.missing = a$num.missing,
        cor.incorrect.m = a$cor.incorrect.m, auc.incorrect.m = a$auc.incorrect.m,
        cor.correct.m = a$cor.correct.m, auc.correct.m = a$auc.correct.m))
  }

  r
}

write.accuracy.tables = function() {
  missing.cells = data.frame(sim.with.various.params())
  write.table(missing.cells, col.names=NA,
    file="git/unmix/missing/sim/missing_cells.txt", sep="\t")

  missing.cells.separate.fractions = data.frame(sim.removing.from.individual.fractions())
  missing.cells.separate.fractions$percent.missing =
    100 * (missing.cells.separate.fractions$num.missing / prod(dim(m1)))
  write.table(missing.cells.separate.fractions, col.names=NA,
    file="git/unmix/missing/sim/missing_cells_separate_fractions.txt", sep="\t")
}

plot.it = function() {
  missing.cells = read.table(
    "git/unmix/missing/sim/missing_cells.txt",
    sep="\t", header=TRUE, row.names=1, as.is=TRUE, stringsAsFactors=FALSE)
  missing.cells.separate.fractions = read.table(
    "git/unmix/missing/sim/missing_cells_separate_fractions.txt",
    sep="\t", header=TRUE, row.names=1, as.is=TRUE, stringsAsFactors=FALSE)

  pdf("git/unmix/missing/sim/accuracy with incorrect matrix.pdf", width=8, height=7)
  par(mfrow=c(2,2))

  plot(missing.cells$num.missing, missing.cells$cor.incorrect.m, ylim=c(0,1), pch=20,
    xlab="number of cells missing", ylab="correlation",
    main="Accuracy with cells missing")

  plot(missing.cells$num.missing, missing.cells$auc.incorrect.m, ylim=c(0.5,1), pch=20,
    xlab="number of cells missing", ylab="area under the curve",
    main="Accuracy with cells missing")


  plot(missing.cells.separate.fractions$percent.missing, missing.cells.separate.fractions$cor.incorrect.m, ylim=c(0,1), pch=20,
    xlab="percent incorrect entries", ylab="correlation",
    main="Accuracy with incorrect sorting")

  plot(missing.cells.separate.fractions$percent.missing, missing.cells.separate.fractions$auc.incorrect.m, ylim=c(0.5,1), pch=20,
    xlab="percent incorrect entries", ylab="area under the curve",
    main="Accuracy with incorrect sorting")

  dev.off()
}

# write.accuracy.tables()
plot.it()

