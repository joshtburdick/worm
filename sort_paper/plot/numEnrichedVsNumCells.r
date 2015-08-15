# Plots the number of genes enriched (or depleted) vs. the number of
# cells in a sorted fraction.

# XXX we assume this was run already
source("git/sort_paper/plot/numEnrichedInFractions.r")

# the sort matrix
source("git/sort_paper/unmix/sortMatrix.r")

n.cells = apply(m.unnormalized >= 0.5, 1, sum)

# number enriched / depleted per fraction
n.enriched = apply(r.sort.only.averaged >= 2, 2, sum)
n.depleted = apply(r.sort.only.averaged <= -2, 2, sum)

s = intersect(names(n.cells), names(n.enriched))

pdf("git/sort_paper/plot/numEnrichedVsNumCells.pdf",
  width=5.5, height=5.5)
par(mar=c(5,5,4,1) + 0.1)
xlim = c(0, max(n.cells[s]))
ylim = c(0, max(c(n.enriched[s], n.depleted[s])))
plot(n.cells[s], n.enriched[s], pch=20, col="#ff0000b0",
  xlim=xlim, ylim=ylim, cex.main=1.1,
  main = "Genes enriched in different\nnumbers of sorted cells",
  xlab = "Number of sorted cells",
  ylab = "Number of genes enriched / depleted")
par(new=TRUE)
plot(n.cells[s], n.depleted[s], pch=20, col="#0000ffb0",
  xlim=xlim, ylim=ylim, main="", xlab="", ylab="")

legend("topright", legend=c("enriched", "depleted"),
  col=c("#ff0000b0", "#0000ffb0"), pch=20)
dev.off()

# tests of significance of trend
cat("enriched Wilcoxon =",
  2 * wilcox.test(n.cells[s], n.enriched[s], paired=TRUE, exact=FALSE)$p.value,
  "\n")

cat("depleted Wilcoxon =",
  2 * wilcox.test(n.cells[s], n.depleted[s], paired=TRUE, exact=FALSE)$p.value,
  "\n")




