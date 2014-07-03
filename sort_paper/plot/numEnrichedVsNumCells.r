# Plots the number of genes enriched (or depleted) vs. the number of
# cells in a sorted fraction.

source("git/sort_paper/plot/numEnrichedInFractions.r")

# the sort matrix
source("git/sort_paper/unmix/sortMatrix.r")

n.cells = apply(m.unnormalized >= 0.5, 1, sum)

# number enriched / depleted per fraction
n.enriched = apply(r.sort.only.averaged >= 2, 2, sum)
n.depleted = apply(r.sort.only.averaged <= -2, 2, sum)


s = intersect(names(n.cells), names(n.enriched))


pdf("git/sort_paper/plot/numEnrichedVsNumCells.pdf")
plot(n.cells[s], n.enriched[s], pch=20, col="#ff000080")
par(new=TRUE)
plot(n.cells[s], n.depleted[s], pch=20, col="#0000ff80")

dev.off()

