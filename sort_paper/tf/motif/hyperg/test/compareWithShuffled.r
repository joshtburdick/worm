# Comparing enrichments with shuffled versions of same.

library("reshape2")

# versions with shuffled motifs
load("git/sort_paper/tf/motif/hyperg/allResults/jolma2013_shuffled/facs_0.892.Rdata")
facs.892.shuffled.motif.enrich = enrich

load("git/sort_paper/tf/motif/hyperg/allResults/jolma2013_shuffled/facs.Rdata")
facs.shuffled.motif.enrich = enrich

load("git/sort_paper/tf/motif/hyperg/allResults/jolma2013_shuffled/hier.300.clusters.Rdata")
hier.300.shuffled.motif.enrich = enrich

m = dimnames(facs.892.shuffled.motif.enrich)[[1]]

# enrichments
load("git/sort_paper/tf/motif/hyperg/allResults/5kb/facs_0.892.Rdata")
facs.892.enrich = enrich[m,,,,,]
load("git/sort_paper/tf/motif/hyperg/allResults/5kb/facs.Rdata")
facs.enrich = enrich[m,,,,,]
load("git/sort_paper/tf/motif/hyperg/allResults/5kb/hier.300.clusters.Rdata")
hier.300.enrich = enrich[m,,,,,]

# useful transformations
# min.p = function(e) apply(e[,,"p",,,], c(1,2), min)
# max.enrich = function(e) apply(e[,,"enrich",,,], c(1,2), max)
p1 = function(e) -log10(e[,,"p",,,])
e1 = function(e) e[,,"enrich",,,]    # log2(1e-6 + e[,,"enrich",,,])


plot.it = function(x, y, main) {
  lim = range(c(x,y))
  plot(x, y, xlim=lim, ylim=lim,
    main=main, xlab="not shuffled", ylab="shuffled",
    pch=183, font=5, cex=2, col="#00000020")
  abline(0,1)
}

if (FALSE) {

  # pdf("git/sort_paper/tf/motif/hyperg/test/compareWithShuffled.pdf",
  #   width=8, height=8)
  png("git/sort_paper/tf/motif/hyperg/test/compareWithShuffled.png",
    width=1500, height=1000)

  par(mfrow=c(2,3))

  plot.it(e1(facs.892.enrich), e1(facs.892.shuffled.motif.enrich),
    main = "Motifs in FACS fractions (cutoff = 0.892), log2 enrichment")
  plot.it(e1(facs.enrich), e1(facs.shuffled.motif.enrich),
    main = "Motifs in FACS fractions (cutoff = 2), log2 enrichment")
  plot.it(e1(hier.300.enrich), e1(hier.300.shuffled.motif.enrich),
    main="300 clusters, log2 enrichment")

  plot.it(p1(facs.892.enrich), p1(facs.892.shuffled.motif.enrich),
    main="Motifs in FACS fractions (cutoff = 0.892), -log10(p)")
  plot.it(p1(facs.enrich), p1(facs.shuffled.motif.enrich),
    main="Motifs in FACS fractions (cutoff = 2), -log10(p)")
  plot.it(p1(hier.300.enrich), p1(hier.300.shuffled.motif.enrich),
    main="300 clusters, -log10(p)")

  dev.off()
}

# table of enrichments before and after clustering
neg.log.p.cluster = melt(hier.300.enrich[,,"p.corr",,,])
neg.log.p.cluster[,6] = -log10(neg.log.p.cluster[,6])
colnames(neg.log.p.cluster)[[6]] = "nlp"
neg.log.p.cluster$shuffled.nlp =
  -log10(melt(hier.300.shuffled.motif.enrich[,,"p.corr",,,])$value)

# motif-cluster pairs which are "more significant before shuffling"
signif.before.shuffling = neg.log.p.cluster[
  neg.log.p.cluster$nlp >= 10 & neg.log.p.cluster$shuffled.nlp <= 1 , ]
print("cases significant before shuffling (but not after)")
print(table(
  paste(signif.before.shuffling[,1], signif.before.shuffling[,2])))

# pairs significant before and after
signif.with.or.without.shuffle = neg.log.p.cluster[
  neg.log.p.cluster$nlp >= 15 & neg.log.p.cluster$shuffled.nlp >= 15 , ]
print("cases significant with or without shuffling")
print(table(
  paste(signif.with.or.without.shuffle[,1],
    signif.with.or.without.shuffle[,2])))


