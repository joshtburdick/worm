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

m1 = m.cell[ rownames(reporters[1:30,]) , ]

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
    m1 = m1 * mask
    num.missing = sum( mask==0 )
  }

  # compute fractions (with the altered matrix)
  sim.frac = expr %*% t(m1)

  # do deconvolving (with incorrect and correct sort matrix)  
  x.d.incorrect.m = sim.frac %*% pseudoinverse(t(m))
  x.d.correct.m = sim.frac %*% pseudoinverse(t(m1))

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

  for(num.lineages.missing in c(1:9))
    for(min.cells in c(1, 20, 40)) {
cat(num.lineages.missing, min.cells, "\n")
      max.cells = min.cells + 20
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

  for(num.lineages.missing in c(1:30))
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

missing.cells = data.frame(sim.with.various.params())
write.table(missing.cells,
  file="git/unmix/missing/sim/missing_cells.txt", sep="\t")
missing.cells.separate.fractions = data.frame(sim.removing.from.individual.fractions())
write.table(missing.cells.separate.fractions,
  file="git/unmix/missing/sim/missing_cells_separate_fractions.txt", sep="\t")

plot.it = function() {

  pdf("git/unmix/missing/sim/missing_cells_1.pdf", width=9, height=3)
  par(mfrow=c(1,3))

  plot(missing.cells$num.missing, missing.cells$cor.incorrect.m, ylim=c(0,1), pch=1,
    xlab="number of cells missing", ylab="correlation",
    main="")
  par(new=TRUE)
  plot(missing.cells$num.missing, missing.cells$cor.correct.m, ylim=c(0,1), pch=20,
    xlab="", ylab="", main="")
  legend("bottomleft", legend = c("correct cell-sorting matrix", "incorrect cell-sorting matrix"),
    pch=c(20, 1))

  plot(missing.cells$num.missing, missing.cells$auc.incorrect.m, ylim=c(0,1), pch=1,
    xlab="number of cells missing", ylab="area under the curve",
    main="")
  par(new=TRUE)
  plot(missing.cells$num.missing, missing.cells$auc.correct.m, ylim=c(0,1), pch=20,
    xlab="", ylab="", main="")
  legend("bottomleft", c("correct cell-sorting matrix", "incorrect cell-sorting matrix"),
    pch=c(20, 1))

  plot(missing.cells$cor.correct.m, missing.cells$cor.incorrect.m,
    xlim=c(0,1), ylim=c(0,1), pch=20, main="",
    xlab="correlation with correct cell-sorting matrix",
    ylab="correlation with incorrect cell-sorting matrix")
  abline(0, 1, col="#303030")

  dev.off()

  pdf("git/unmix/missing/sim/missing_cells_2.pdf", width=9, height=3)
  par(mfrow=c(1,3))

  plot(missing.cells.separate.fractions$num.missing, missing.cells.separate.fractions$cor.incorrect.m, ylim=c(0,1), pch=1,
    xlab="number of cells missing", ylab="correlation",
    main="")
  par(new=TRUE)
  plot(missing.cells.separate.fractions$num.missing, missing.cells.separate.fractions$cor.correct.m, ylim=c(0,1), pch=20,
    xlab="", ylab="", main="")
  legend("bottomleft", legend = c("correct cell-sorting matrix", "incorrect cell-sorting matrix"),
    pch=c(20, 1))

  plot(missing.cells.separate.fractions$num.missing, missing.cells.separate.fractions$auc.incorrect.m, ylim=c(0,1), pch=1,
    xlab="number of cells missing", ylab="area under the curve",
    main="")
  par(new=TRUE)
  plot(missing.cells.separate.fractions$num.missing, missing.cells.separate.fractions$auc.correct.m, ylim=c(0,1), pch=20,
    xlab="", ylab="", main="")
  legend("bottomleft", c("correct cell-sorting matrix", "incorrect cell-sorting matrix"),
    pch=c(20, 1))

  plot(missing.cells.separate.fractions$cor.correct.m, missing.cells.separate.fractions$cor.incorrect.m,
    xlim=c(0,1), ylim=c(0,1), pch=20, main="",
    xlab="correlation with correct cell-sorting matrix",
    ylab="correlation with incorrect cell-sorting matrix")
  abline(0, 1, col="#303030")

  dev.off()
}

plot.it()

