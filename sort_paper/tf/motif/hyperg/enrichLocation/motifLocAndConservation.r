# Plots how much motifs were enriched at particular
# locations, or with higher conservation.

source("git/utils.r")

source("git/sort_paper/tf/motif/hughes/motifInfo.r")

# for filtering motifs
load("git/sort_paper/tf/motif/hughes/motifCluster.Rdata")
motif.clustering = cutree(hughes.motif.cluster[["all"]], h=0.01)

motif.ortholog.summary = sapply(orthologs.by.motif,
  function(a) {
    if (length(a) <= 21)
      paste(a, collapse=" ")
    else
      paste(paste(a[1:20], collapse=" "),
        "and", length(a) - 20, "others", collapse=" ")
  })

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

a = a[ order(a$lp.corr, decreasing=TRUE) , ]

# filter redundant motifs
a = a[ !duplicated(paste(motif.clustering[ a$motif ], a$group)) , ]

# add on motif annotation
a$motif.name = motif.info[ a$motif, "motif.name" ]
a$gene.or.ortholog = motif.info[ a$motif, "gene" ]
m = a$motif %in% names(motif.ortholog.summary)
a[m, "gene.or.ortholog"] = motif.ortholog.summary[ a[m,"motif"] ]
a$gene.or.ortholog[ is.na(a$gene.or.ortholog) ] = ""

# write out data used for plotting
write.tsv(a,
  "git/sort_paper/tf/motif/hyperg/enrichLocation/motifLocAndConservation.tsv")

# Plots each of these four categories separately, compared
# to motif enrichment.
plot.1 = function() {
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
}

# Plots each of these together, colored by motif enrichment.
plot.2 = function(coloring) {
  a = a[ !is.na(a$lp.cons) , ]
  a = a[ is.finite(a$lp.dist) , ]
  a$lp.corr.1 = pmin(a$lp.corr / 10, 1)
  a$lp.dist.1 = ifelse(a$dist.closer, a$lp.dist, -a$lp.dist)
  a$lp.cons.1 = ifelse(a$cons.higher, a$lp.cons, -a$lp.cons)

  plot(a$lp.dist.1, a$lp.cons.1, xlim=c(-305, 305),
    ylim=c(-40,40),
    xlab="Distance to TSS", ylab="Conservation",
    xaxt="n", yaxt="n",
    pch=183, font=5, cex=0.7, col=coloring(a$lp.corr.1))
  abline(h=0, col="#00000080")
  abline(v=0, col="#00000080")
  a = 100*c(-3:3)
  axis(1, at=a, labels=abs(a), cex=0.8)
  a = 20*c(-2:2)
  axis(2, at=a, labels=abs(a), cex=0.8)

  legend("topleft", title="Enrichment significance", pch=20,
#    legend=c(4,6,8,10),
    legend=c(expression(10^{-4}),
      expression(10^{-6}),
      expression(10^{-8}),
      expression("" < 10^{-10})),
      col=coloring(c(2:5)/5), cex=0.7)
    mtext("Further", side=1, adj=0, line=2)
    mtext("Closer", side=1, adj=1, line=2)
    mtext("Lower", side=2, adj=0, line=2)
    mtext("Higher", side=2, adj=1, line=2)
}

pdf("git/sort_paper/tf/motif/hyperg/enrichLocation/motifLocAndConservation_BW.pdf",
  width=10, height=7.5)
plot.2(function(x) hsv(0, 0, 1-x, 0.8))
dev.off()

pdf("git/sort_paper/tf/motif/hyperg/enrichLocation/motifLocAndConservation_C.pdf",
  width=10, height=7.5)
plot.2(function(x) hsv(1-x, 1, 1, 0.8))
dev.off()


