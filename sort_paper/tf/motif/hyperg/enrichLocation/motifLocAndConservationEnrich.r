# Plots how much motifs were enriched at particular
# locations, or with higher conservation, using enrichment
# (based on a cutoff.)

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

r = read.tsv(gzfile("git/cluster/motif/plot/distConservationEnrich2.tsv.gz"))
# r = r[1:1000,]    # not finished yet

# possibly only keep presumably significant enrichments
significant = (r$cons.p.corr <= 0.05) & (r$dist.p.corr <= 0.05) 
cat("proportion with cons and dist significant =", mean(significant), "\n")
# r = r[ significant , ]

e1 = enrich[ , c("motif", "group", "enrich", "p.corr") ]
r1 = r[ , c("motif", "group", "dist.enrich", "cons.enrich", "dist.p.corr", "cons.p.corr") ]

a = merge(e1, r1)
a = a[ a$p.corr <= 0.05 , ]
a[,3] = log2(a[,"enrich"])
colnames(a)[3] = "log.enrich"
a[,4] = -log10(a[,"p.corr"])
colnames(a)[4] = "lp.corr"

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
  "git/sort_paper/tf/motif/hyperg/enrichLocation/motifLocAndConservationWithEffectSize.tsv")

# plot most significant last
a = a[ nrow(a):1 , ]

# Plots each of these together, colored by motif enrichment p-value.
plot.2 = function(coloring) {
#  lim = c(0.3, 0.7)

  a$log.enrich.1 = pmin(a$log.enrich, 1) 
  a$lp.corr.1 = pmin(a$lp.corr / 10, 1)
#  a$dist.effect.size = pmin(pmax(a$dist.effect.size, lim[1]), lim[2])
#  a$cons.effect.size = pmin(pmax(a$cons.effect.size, lim[1]), lim[2])
# browser()
  plot(a$dist.enrich, a$cons.enrich,
    xlab="Enrichment of motifs within 1kb of 5' end of transcript",
    ylab="Enrichment of motifs with conservation > 0.5",
    xaxt="n", yaxt="n",
    pch=183, font=5, cex=0.7, col=coloring(a$lp.corr.1))
  abline(h=1, col="#00000080")
  abline(v=1, col="#00000080")
  a = c(3:7) / 10
  axis(1, cex=0.8)   # at=a, labels=abs(a))
  axis(2, cex=0.8)   # at=a, labels=abs(a))

if (TRUE) {
  legend("bottomright", title="Motif enrichment", pch=20,
    legend=c(expression(10^{-4}),
      expression(10^{-6}),
      expression(10^{-8}),
      expression("" < 10^{-10})),
    col=coloring(c(2:5)/5), cex=0.7)
}
}

# Plots each of these together, colored by motif enrichment p-value.
plot.colored.by.enrichment = function() {

  a$log.enrich.1 = pmin(a$log.enrich, 1) 
  a$lp.corr.1 = pmin(a$lp.corr / 10, 1)

  # test of whether either of these is significant,
  # for shading the points
  sig = (a$dist.p.corr <= 0.05) | (a$cons.p.corr <= 0.05)

  # how to color dots
  coloring = function(x, e)
    hsv(1-x, ifelse(e, 1, 0.5), 0.8, 0.8)

  plot(a$dist.enrich, a$cons.enrich,   # xlim=lim, ylim=lim,
    xlab="Enrichment of motifs within 1kb of 5' end of transcript",
    ylab="Enrichment of motifs with conservation > 0.5",
    xaxt="n", yaxt="n",
    pch=183, font=5, cex=0.7,
#    col=coloring(a$log.enrich.1))
    col=coloring(a$log.enrich.1, sig))
  abline(h=1, col="#00000080")
  abline(v=1, col="#00000080")
  a = c(3:7) / 10
  axis(1, cex=0.8)   # at=a, labels=abs(a))
  axis(2, cex=0.8)   # at=a, labels=abs(a))

if (TRUE) {
  legend("bottomleft", ncol=2,
    title=expression(log[2] * "(motif enrichment near genes in cluster)"),
    pch=20, legend=c(c(10,8,6,4,2) / 10, paste(c(10,8,6,4,2)/10, "*")),
    col=c(coloring(c(5:1)/5, FALSE), coloring(c(5:1)/5, TRUE)), cex=0.62)
#  legend("bottomleft",
#    title="", pch=20, bty="n", legend = rep("", 4),
#    col=coloring(c(5:2)/5, TRUE), cex=0.68)
}

}

# Same, but somewhat simplified. This doesn't grey out
# non-significant things.
plot.colored.by.enrichment.simplified = function() {

  a$log.enrich.1 = pmin(a$log.enrich, 1) 
  a$lp.corr.1 = pmin(a$lp.corr / 10, 1)

  # test of whether either of these is significant,
  # for shading the points
  sig = (a$dist.p.corr <= 0.05) | (a$cons.p.corr <= 0.05)

  # how to color dots (ignoring shading, for now)
  coloring = function(x, e)
    hsv(1-x, 1, 0.8, 0.8)

  plot(a$dist.enrich, a$cons.enrich,   # xlim=lim, ylim=lim,
    xlab="Enrichment of motifs within 1kb of 5' end of transcript",
    ylab="Enrichment of motifs with conservation > 0.5",
    xaxt="n", yaxt="n",
    pch=183, font=5, cex=0.7,
#    col=coloring(a$log.enrich.1))
    col=coloring(a$log.enrich.1, sig))
  abline(h=1, col="#00000080")
  abline(v=1, col="#00000080")
  a = c(3:7) / 10
  axis(1, cex=0.8)   # at=a, labels=abs(a))
  axis(2, cex=0.8)   # at=a, labels=abs(a))

if (TRUE) {
  legend("bottomleft",
    title=expression(log[2] * "(motif enrichment)"),
    pch=20, legend=c(c(10,8,6,4,2) / 10),
    col=coloring(c(5:1)/5, TRUE), cex=0.87)
#  legend("bottomleft",
#    title="", pch=20, bty="n", legend = rep("", 4),
#    col=coloring(c(5:2)/5, TRUE), cex=0.68)
  a = 0.9
    mtext("Further from 5' end", side=1, adj=0, line=2, cex=a)
    mtext("Closer to 5' end", side=1, adj=1, line=2, cex=a)
    mtext("Less conserved", side=2, adj=0, line=2, cex=a)
    mtext("More conserved", side=2, adj=1, line=2, cex=a)
}

}



if (FALSE) {
  pdf("git/sort_paper/tf/motif/hyperg/enrichLocation/motifLocAndConservationEnrich colored by p.pdf",
    width=7.5, height=7.5)
  plot.2(function(x) hsv(1-x, 1, 0.8, 0.5))
  dev.off()

pdf("git/sort_paper/tf/motif/hyperg/enrichLocation/motifLocAndConservationEnrich.pdf",
  width=7.5, height=7.5)
plot.colored.by.enrichment()
dev.off()
}

pdf("git/sort_paper/tf/motif/hyperg/enrichLocation/motifLocAndConservationEnrichSimplified.pdf",
  width=7.5, height=7.5)
plot.colored.by.enrichment.simplified()
dev.off()


