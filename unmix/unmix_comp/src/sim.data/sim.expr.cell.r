# Simulated lineage-based patterns.

load("../data/tree_utils.Rdata")
source("sim.data/sim.expr.r")

# generate simulated data with one lineage on
cell.lineage.1 = cell.lineage.matrix[apply(cell.lineage.matrix, 1, sum)>=5,]
set.seed(42)
sim.expr.cell = sim.expr(cell.lineage.1)

save(sim.expr.cell, file="sim.data/sim.expr.cell.Rdata")

