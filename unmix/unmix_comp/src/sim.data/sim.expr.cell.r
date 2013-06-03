# Simulated lineage-based patterns.

load("../data/tree_utils.Rdata")
source("sim.data/sim.expr.r")
source("sim.data/sim.expr.gamma.r")

# generate simulated data with one lineage on
cell.lineage.1 = cell.lineage.matrix[apply(cell.lineage.matrix, 1, sum)>=5,]
set.seed(42)
sim.expr.cell = sim.expr(cell.lineage.1)

save(sim.expr.cell, file="sim.data/sim.expr.cell.Rdata")

# gamma-distributed version of this
sim.expr.cell.gamma = sim.expr.gamma(cell.lineage.1)
save(sim.expr.cell.gamma, file="sim.data/sim.expr.cell.gamma.Rdata")

