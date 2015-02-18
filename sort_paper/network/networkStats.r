# Computes various network stats.

if (TRUE) {
load("git/sort_paper/tf/motif/hyperg/allResults/5kb/hier.300.clusters.Rdata")
motif1.min.p = apply(enrich[,,"p.corr",,,], c(1,2), min)

load("git/sort_paper/tf/motif/hyperg/allResults/hughes/hier.300.clusters.Rdata")
hughes.min.p = apply(enrich[,,"p.corr",,,], c(1,2), min) 

motif.min.p = rbind(hughes.min.p, motif1.min.p)

load("git/sort_paper/tf/motif/hyperg/allResults/chip/hier.300.clusters.Rdata")
chip.min.p = apply(enrich[,,"p.corr",,,], c(1,2), min) 
}

pdf("git/sort_paper/network/networkStats.pdf",
  width=10, height=6)
par(mfrow=c(2,4))

hist(apply(motif.min.p <= 0.05, 1, sum),
  main="Clusters per motif",
  xlab="Number of clusters",
  breaks=50, col="red")
mtext("p <= 0.05", 3, cex=0.7)
hist(apply(motif.min.p <= 1e-5, 1, sum),
  main="Clusters per motif",
  xlab="Number of clusters",
  breaks=50, col="red")
mtext("p <= 1e-5", 3, cex=0.7)
hist(apply(chip.min.p <= 0.05, 1, sum),
  main="Clusters per ChIP signal",
  xlab="Number of clusters",
  breaks=25, col="red")
mtext("p <= 0.05", 3, cex=0.7)
hist(apply(chip.min.p <= 1e-5, 1, sum),
  main="Clusters per ChIP signal",
  xlab="Number of clusters",
  breaks=25, col="red")
mtext("p <= 1e-5", 3, cex=0.7)

hist(apply(motif.min.p <= 0.05, 2, sum),
  main="Motifs per cluster",
  xlab="Number of motifs",
  breaks=50, col="blue")
mtext("p <= 0.05", 3, cex=0.7)
hist(apply(motif.min.p <= 1e-5, 2, sum),
  main="Motifs per cluster",
  xlab="Number of motifs",
  breaks=50, col="blue")
mtext("p <= 1e-5", 3, cex=0.7)
hist(apply(chip.min.p <= 0.05, 2, sum),
  main="ChIP signals per cluster",
  xlab="Number of ChIP signals",
  breaks=25, col="blue")
mtext("p <= 0.05", 3, cex=0.7)
hist(apply(chip.min.p <= 1e-5, 2, sum),
  main="ChIP signals per cluster",
  xlab="Number of ChIP signals",
  breaks=25, col="blue")
mtext("p <= 1e-5", 3, cex=0.7)

dev.off()

