# Computes "cell-sorting matrix", estimating how much each cell
# was present in each sorted fraction.

# cell volume (and weight)
source("git/unmix/image/cell_volume.r")

# FIXME this probably could be improved (also, it should be renamed)
gene.on.off = read.csv("R/unmix/sort_paper/unmix/image/sort.matrix.csv.gz",
  as.is=TRUE, row.names=1)

sort.matrix = t( t(gene.on.off) * cell.weights[colnames(gene.on.off),"w"] )
sort.matrix[,"P0"] = 0

write.table(sort.matrix, file=gzfile("git/unmix/image/sort.matrix.tsv.gz"),
  sep="\t")

save(sort.matrix, file="git/unmix/image/sort_matrix.Rdata")

