# Computes "cell-sorting matrix", estimating how much each cell
# was present in each sorted fraction.

# cell volume (and weight)
source("git/unmix/image/cell_volume.r")

# FIXME this probably could be improved (also, it should be renamed)
gene.on.off = read.csv("R/unmix/sort_paper/unmix/image/sort.matrix.csv.gz",
  as.is=TRUE, row.names=1)

# compute "negative" portions
gene.on.off = {
  negative.fractions = 1 - gene.on.off
  rownames(negative.fractions) = paste(rownames(gene.on.off), "_minus", sep="")
  rbind(gene.on.off, negative.fractions[-1,])
}

# version of this, not weighted by cell volume
sort.matrix.unweighted = as.matrix(gene.on.off)
sort.matrix.unweighted[,"P0"] = 0

# same, weighted by cell volume
sort.matrix = t( t(gene.on.off) * cell.weights[colnames(gene.on.off),"w"] )
sort.matrix[,"P0"] = 0

write.table(sort.matrix, file=gzfile("git/unmix/image/sort.matrix.tsv.gz"),
  sep="\t")

save(sort.matrix, sort.matrix.unweighted, file="git/unmix/image/sort_matrix.Rdata")

