# Counts motifs enriched upstream of each cluster,
# for different clustering methods.

source("git/utils.r")

cluster.dir = "git/cluster/hierarchical/"

# Gets statistics for chi-squared version of this.
get.chisquared.counts = function(f, filter.enrich=FALSE) {
  r = read.table(f, sep="\t", header=TRUE, row.names=NULL, as.is=TRUE)
  if (filter.enrich)
    r = r[ r$enrichment > 1 , ]
  r = r[ r$p.fdr <= 0.05 , ]

  data.frame(method="chi.sq", num.pairs = nrow(r),
    num.clusters = length(unique(r$cluster)),
    num.motifs = length(unique(r$motif)))
}

# Gets statistics for the Mann-Whitney Wilcox version of this test.
get.wilcox.counts = function(f, filter.enrich=FALSE) {
  r = read.table(f, sep="\t", header=TRUE, row.names=NULL, as.is=TRUE)
  if (filter.enrich)
    r = r[ r$motifs.cluster.mean > r$motifs.background.mean , ]
  r = r[ r$p.fdr <= 0.05 , ]

  data.frame(method="wilcox", num.pairs = nrow(r),
    num.clusters = length(unique(r$cluster)),
    num.motifs = length(unique(r$motif)))
}

num.motifs.enriched = NULL
for(d in list.files(cluster.dir)[c(1:5,7:11)]) {
  cat(paste(backspace.string, d))

  # first, including the actual motifs
  chisq = get.chisquared.counts(paste(cluster.dir,
    d, "/uniqueKnownMotifEnrichment_5kb.tsv", sep=""))
  chisq.cons = get.chisquared.counts(paste(cluster.dir,
    d, "/uniqueKnownMotifEnrichment_5kb_0.5cons.tsv", sep=""))
  wilcox = get.wilcox.counts(paste(cluster.dir,
    d, "/WilcoxUniqueKnownMotifEnrichment_5kb.tsv",sep=""))
  wilcox.cons = get.wilcox.counts(paste(cluster.dir,
    d, "/WilcoxUniqueKnownMotifEnrichment_5kb_0.5cons.tsv", sep=""))

  r1 = rbind(
    cbind(data.frame(method="chi.sq", conservation="no"), chisq),
    cbind(data.frame(method="chi.sq", conservation="yes"), chisq.cons),
    cbind(data.frame(method="wilcox", conservation="no"), wilcox),
    cbind(data.frame(method="wilcox", conservation="yes"), wilcox.cons))

  num.motifs.enriched =
    rbind(num.motifs.enriched, cbind(clustering=d, shuffled="no", r1))

  # also tack on random motifs
  chisq = get.chisquared.counts(paste(cluster.dir,
    d, "/uniqueShuffledKnownMotifEnrichment_5kb.tsv", sep=""))
  chisq.cons = get.chisquared.counts(paste(cluster.dir,
    d, "/uniqueShuffledKnownMotifEnrichment_5kb_0.5cons.tsv", sep=""))
  wilcox = get.wilcox.counts(paste(cluster.dir,
    d, "/WilcoxUniqueShuffledKnownMotifEnrichment_5kb.tsv",sep=""))
  wilcox.cons = get.wilcox.counts(paste(cluster.dir,
    d, "/WilcoxUniqueShuffledKnownMotifEnrichment_5kb_0.5cons.tsv", sep=""))

  r1 = rbind(
    cbind(data.frame(method="chi.sq", conservation="no"), chisq),
    cbind(data.frame(method="chi.sq", conservation="yes"), chisq.cons),
    cbind(data.frame(method="wilcox", conservation="no"), wilcox),
    cbind(data.frame(method="wilcox", conservation="yes"), wilcox.cons))

  num.motifs.enriched =
    rbind(num.motifs.enriched, cbind(clustering=d, shuffled="yes", r1))


}

cat("\n")

num.motifs.enriched =
  num.motifs.enriched[ order(num.motifs.enriched$method) , ]

write.tsv(num.motifs.enriched, "git/cluster/motif/num.motifs.enriched.tsv")

