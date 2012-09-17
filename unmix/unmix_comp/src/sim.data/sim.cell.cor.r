# Simulated data, based on the correlation matrix.

library("corpcor")
library("mnormt")

# calculate covariance of movies
expr.cell = as.matrix(read.table("../data/exprCell.tsv",
  sep="\t", header=TRUE, row.names=1, as.is=TRUE))
cor.cell.shrunk = cor.shrink(expr.cell, lambda=0.05)

# simulated data
set.seed(42)
sim.cell.cor = rmnorm(200, varcov=cor.cell.shrunk)
sim.cell.cor[sim.cell.cor<0] = 0
rownames(sim.cell.cor) = paste("corr. sim.", 1:dim(sim.cell.cor)[1])

save(sim.cell.cor, file="sim.data/sim.cell.cor.Rdata")

