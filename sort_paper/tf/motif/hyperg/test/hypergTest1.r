# More trying to understand why everything appears significant.

load("git/sort_paper/tf/motif/hyperg/allResults/5kb/hier.300.clusters.Rdata")

p = enrich[,,"p",,,]
p.min = apply(enrich[,,"p",,,], c(1,2), min)

e = enrich[,,"enrich",,,]
e[is.na(e)] = 0



p.sig = 1 * (p <= 0.05)
g.cluster = apply(enrich[,,"g.cluster",,,], c(1,2), max)
m.total = apply(enrich[,,"m.total",,,], c(1,2), max)

e.max = apply(e,c(1,2), max)

if (TRUE) {
pdf("git/sort_paper/tf/motif/hyperg/test/hypergTest1.pdf",
  width=7.5, height=10)
par(mfrow=c(2,2))

  hist(log2(e), breaks=100, col="grey",
    xlab="log2(motif enrichment)")
  plot.new()
  hist(log2(apply(e,c(1,2), min)), breaks=100, col="grey",
    xlab="minimum log2(motif enrichment)")
  hist(log2(apply(e,c(1,2), max)), breaks=100, col="grey",
    xlab="maximum log2(motif enrichment)")

dev.off()
}


# 
if (FALSE) {
plot(g.cluster, -log10(p.min), pch=20, col="#00000080")
plot(m.total, -log10(p.min), pch=20, col="#00000080")

plot(g.cluster, e.max, pch=20, col="#00000080")
plot(m.total, e.max, pch=20, col="#00000080")
}

