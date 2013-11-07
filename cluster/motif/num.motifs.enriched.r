# Counts motifs enriched upstream of each cluster,
# for different clustering methods.

source("git/utils.r")

cluster.dir = "git/cluster/hierarchical/"

# Gets statistics for the Mann-Whitney U version of this test.
get.mw.counts = function(f, filter.enrich) {
#  f = paste(f, "/MW_knownMotifEnrichment_5kb.tsv", sep="")
# cat("\n", f, "\n")
  r = read.tsv(f)
  if (filter.enrich)
    r = r[ r$a.mean > r$b.mean , ]
  r = r[ r$p.bh <= 0.05 , ]

  list(num.pairs = nrow(r),
    num.clusters = length(unique(r$cluster)),
    num.motifs = length(unique(r$name)))
}

# Gets statistics for chi-squared version of this.
get.chisquared.counts = function(f, filter.enrich) {
#  f = paste(f, "/knownMotifEnrichment_5kb.tsv", sep="")
# cat("\n", f, "\n")
  r = read.tsv(f)
  if (filter.enrich)
    r = r[ r$enrichment > 1 , ]
  r = r[ r$p.fdr <= 0.05 , ]

  list(num.pairs = nrow(r),
    num.clusters = length(unique(r$cluster)),
    num.motifs = length(unique(r$motif)))
}

num.motifs.enriched = NULL
for(d in list.files(cluster.dir)) {
  cat(paste(backspace.string, d))
  mw = get.mw.counts(paste(cluster.dir,
    d, "/MW_knownMotifEnrichment_5kb.tsv",sep=""))
#  mw.cons = get.mw.counts(paste(cluster.dir,
#    d, "/MW_knownMotifEnrichment_5kb_0.5cons.tsv", sep=""))
  chisq = get.chisquared.counts(paste(cluster.dir,
    d, "/knownMotifEnrichment_5kb.tsv", sep=""))
#  chisq.cons = get.chisquared.counts(paste(cluster.dir,
#    d, "/knownMotifEnrichment_5kb_0.5cons.tsv", sep=""))

  num.motifs.enriched = rbind(num.motifs.enriched,
    data.frame(clustering=d, method="MW", num.pairs=mw$num.pairs,
      num.clusters=mw$num.clusters, num.motifs=mw$num.motifs),
    data.frame(clustering=d, method="chisq", num.pairs=chisq$num.pairs,
      num.clusters=chisq$num.clusters, num.motifs=chisq$num.motifs))
#    data.frame(clustering=d, method="MW.cons", num.pairs=mw.cons$num.pairs,
#      num.clusters=mw.cons$num.clusters, num.motifs=mw.cons$num.motifs))
#    data.frame(clustering=d, method="chisq.cons", num.pairs=chisq.cons$num.pairs,
#      num.clusters=chisq.cons$num.clusters, num.motifs=chisq.cons$num.motifs))
}
cat("\n")
write.tsv(num.motifs.enriched, "git/cluster/motif/num.motifs.enriched.tsv")

