# Optimizes the pseudocount, and enrichment cutoff,
# for considering genes "enriched" (or "depleted".)

# note: this re-generates the "readRatios.tsv" table
source("git/cluster/readRatios.r")

# Computes how many genes were reproducibly enriched,
# using different pseudocounts and cutoffs.
# Args:
#   r - the "reads per million" counts
#   samples - names of samples (without "(+)" or "(-)")
#   pseudocount - the pseudocount of reads per million to add
#   cutoff - cutoff to use for things being enriched/depleted
compute.enriched = function(r, samples, pseudocount, cutoff) {
  r.pos = r[ , paste(samples, "(+)") ]
  r.neg = r[ , paste(samples, "(-)") ]

  enrich = log2( pseudocount + r.pos ) - log2( pseudocount + r.neg )  

  c( enriched = table( apply(enrich >= cutoff, 1, sum) )[2:4],
    depleted = table( apply(enrich <= -cutoff, 1, sum) )[2:4])
}

# Computes above counts for various cutoffs.
enrichment.multiple.settings = function(r, samples) {
  a = NULL

  for(pseudocount in c(1:10))
    for(cutoff in c(1:6) / 2) {
      cat(backspace.string, pseudocount, cutoff)
      a1 = compute.enriched(r, samples, pseudocount, cutoff)
      a = rbind(a, c(pseudocount=pseudocount, cutoff=cutoff, a1))
    }

  a[ is.na(a) ] = 0
  a
}

cnd1.samples = c("cnd-1 12/14", "cnd-1 1/4", "cnd-1 8/19")
pha4.samples = c("pha-4 12/9", "pha-4 5/9", "pha-4 9/1")


foo = enrichment.multiple.settings(readsPerMillion, cnd1.samples)

es = rbind(
  data.frame(gene="cnd-1", enrichment.multiple.settings(readsPerMillion, cnd1.samples)),
  data.frame(gene="pha-4", enrichment.multiple.settings(readsPerMillion, pha4.samples)))

es$enriched.agreement = round(es$enriched.3 /
  (es$enriched.1 + es$enriched.2 + es$enriched.3), 3)
es$depleted.agreement = round(es$depleted.3 /
  (es$depleted.1 + es$depleted.2 + es$depleted.3), 3)
es$avg.agreement = (es$enriched.agreement + es$depleted.agreement) / 2
es = es[order(es$enriched.agreement, decreasing=TRUE),]

write.tsv(es, "git/sort_paper/enrichment/pseudocountAndCutoffOptimize.tsv")


