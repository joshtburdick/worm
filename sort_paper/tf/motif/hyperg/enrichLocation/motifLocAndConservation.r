# Plots how much motifs were enriched at particular
# locations, or with higher conservation.

source("git/utils.r")

enrich = read.tsv(gzfile("git/sort_paper/tf/motif/hyperg/summary/hughes/hier.300.clusters.tsv.gz"))

r = read.tsv(gzfile("git/cluster/motif/plot/distConservationWilcoxon.tsv.gz"))

e1 = enrich[ , c("motif", "group", "p.corr") ]
r1 = r[ , c("motif", "group", "dist.closer", "dist.p.corr",
  "cons.higher", "cons.p.corr") ]

a = merge(e1, r1)
a[,3] = -log10(a$p.corr)
colnames(a)[3] = "lp.corr"
a[,5] = -log10(a$dist.p.corr)
colnames(a)[5] = "lp.dist"
a[,7] = -log10(a$cons.p.corr)
colnames(a)[7] = "lp.cons"

png("git/sort_paper/tf/motif/hyperg/enrichLocation/motifLocAndConservation.png",
  width=1200, height=600)
par(mfrow=c(2,4))

hist(a[a$dist.closer, "lp.dist"], main="Closer to TSS",
  xlab="Location -log10(p)", breaks=100, col="grey")
hist(a[!a$dist.closer, "lp.dist"], main="Further from TSS",
  xlab="Location -log10(p)", breaks=100, col="grey")
hist(a[a$cons.higher, "lp.cons"], main="More conserved",
  xlab="Conservation -log10(p)", breaks=100, col="grey", xlim=c(0,40))
hist(a[!a$cons.higher, "lp.cons"], main="Less conserved",
  xlab="Conservation -log10(p)", breaks=100, col="grey", xlim=c(0,40))

i = a$dist.closer
plot(a[i,"lp.dist"], a[i,"lp.corr"], pch=20, col="#00000040",
  xlab="Location -log10(p)", ylab="Enrichment -log10(p)",
  main="Closer to TSS", ylim=c(0,70))
plot(a[!i,"lp.dist"], a[!i,"lp.corr"], pch=20, col="#00000040",
  xlab="Location -log10(p)", ylab="Enrichment -log10(p)",
  main="Further from TSS", ylim=c(0,70))

i = a$cons.higher
plot(a[i,"lp.cons"], a[i,"lp.corr"], pch=20, col="#00000040",
  xlab="Conservation -log10(p)", ylab="Enrichment -log10(p)",
  main="More conserved", xlim=c(0,40), ylim=c(0,70))
plot(a[!i,"lp.cons"], a[!i,"lp.corr"], pch=20, col="#00000040",
  xlab="Conservation -log10(p)", ylab="Enrichment -log10(p)",
  main="Less conserved", xlim=c(0,40), ylim=c(0,70))

dev.off()

