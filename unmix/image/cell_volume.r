# Estimates cell volume, and proportion of embryos, just based on the lineage.

# ??? this may not be the lineage to use
cell.time.on.off = read.csv("data/image/cellTimeOnOff.csv.gz",
  header=TRUE, row.names=1, as.is=TRUE)

# I'm not sure why this isn't showing up as 550
cell.time.on.off["MSpppaaa","off"] = 550

# arbitrarily adding a tiny bit of time to very late cells
# this avoids zeros in the weight matrix (except for P0)
# It should only affect ABp[l/r]apapaa[a/p]
cell.time.on.off[ cell.time.on.off$on==550 & cell.time.on.off$off==550, "off"] = 551

load("R/lineage/tree_utils.Rdata")

# depth of a cell, with P0 counted as depth 0, P1 and AB as depth 1, etc.
lineage.depth = apply(cell.lineage.matrix, 2, sum) - 1

cell.weights = cell.time.on.off

# estimated volume of each cell
cell.weights$volume = 2 ^ (-lineage.depth[rownames(cell.weights)])

# weight, as volume * time
cell.weights$w = (1 + cell.weights$off - cell.weights$on) * cell.weights$volume
cell.weights$w = cell.weights$w / sum(cell.weights$w)




