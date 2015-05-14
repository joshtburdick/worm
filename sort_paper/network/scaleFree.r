# Tests if transcriptional network is scale-free, at least
# superficially.

if (TRUE) {
load("git/sort_paper/tf/motif/hyperg/allResults/5kb/hier.300.clusters.Rdata")
motif1.min.p = apply(enrich[,,"p.corr",,,], c(1,2), min)

load("git/sort_paper/tf/motif/hyperg/allResults/hughes/hier.300.clusters.Rdata")
hughes.min.p = apply(enrich[,,"p.corr",,,], c(1,2), min) 

motif.min.p = rbind(hughes.min.p, motif1.min.p)

load("git/sort_paper/tf/motif/hyperg/allResults/chip/hier.300.clusters.Rdata")
chip.min.p = apply(enrich[,,"p.corr",,,], c(1,2), min) 
}



scale.free.regress = function(x, main) {

  h = hist(apply(x, 1, sum), breaks=15, col="grey", main=main, xlab="degree")

  ld = log(h$breaks[-1])
  lc = log(h$counts+1)
  plot(ld, lc, pch=20, xlab="log(degree distribution)", ylab="log(count+1)")
  r = lm(lc ~ ld)
  abline(r$coefficients, lwd=2, col="#aa000080")
  mtext(paste("p = ", signif(summary(r)$coef[2,4], 2)), cex=0.7)

  list(h = h, r = r)  
}


pdf("git/sort_paper/network/scaleFree.pdf", width=10, height=5)
par(mfcol=c(2,4))
if (FALSE) {        # skipping these for now
for (min.p in 10^c(-1:-7)) {
# cutoffs were 1e-3, 1e-3, 1e-4, and 1e-2, respectively
  scale.free.regress(motif.min.p <= min.p, "Motifs")
  scale.free.regress(t(motif.min.p <= min.p), "Clusters (motifs)")
  scale.free.regress(chip.min.p <= min.p, "ChIP signals")
  scale.free.regress(t(chip.min.p <= min.p), "Clusters (ChIP signals)")
}
}
# plot just the most significant
scale.free.regress(motif.min.p <= 1e-4, "Motifs")
scale.free.regress(t(motif.min.p <= 1e-4), "Clusters (motifs)")
scale.free.regress(chip.min.p <= 1e-7, "ChIP signals")
scale.free.regress(t(chip.min.p <= 1e-7), "Clusters (ChIP signals)")

dev.off()

