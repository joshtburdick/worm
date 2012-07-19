# Estimates cell volume, and proportion of embryos, just based on the lineage.

# ??? use the standardized lineage?
cell.time.on.off = read.csv("data/image/cellTimeOnOff.csv.gz",
  header=TRUE, row.names=1, as.is=TRUE)

load("R/lineage/tree_utils.Rdata")

# depth of a cell, with P0 counted as depth 0, P1 and AB as depth 1, etc.
lineage.depth = apply(cell.lineage.matrix, 2, sum) - 1

cell.weights = cell.time.on.off


# cutoff time for including cells
t.cutoff = 580
cell.weights = cell.weights[ cell.weights$on < t.cutoff , ]
cell.weights[ cell.weights$off > t.cutoff , "off" ] = t.cutoff

# estimated volume of each cell
cell.weights$volume = 2 ^ (-lineage.depth[rownames(cell.weights)])

# weight, as volume * time
cell.weights$w = (cell.weights$off - cell.weights$on) * cell.weights$volume
cell.weights$w = cell.weights$w / sum(cell.weights$w)




