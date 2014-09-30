# Computes sizes of upstream regions at different conservation cutoffs.

source("git/utils.r")

# number of upstream bp with different levels of conservation
upstream.cons.dist = list()
for(i in 1:5)
  upstream.cons.dist[[i]] =
    read.tsv(paste0("git/tf/motif/conservation/cons_hist_WS220_",
      i, "kb_upstream.tsv.gz"))

genes = rownames(upstream.cons.dist[[1]])

upstream.dist.kb.levels = 1:5
conservation.levels = c(0, 0.5, 0.7, 0.9)

# cache of upstream conservation at different cutoffs,
# indexed first by upstream distance, then conservation
upstream.bp.at.cutoff = array.from.dimnames(list(
  "upstream.dist.kb" = upstream.dist.kb.levels,
  "conservation" = conservation.levels,
  "gene" = genes))

for(upstream.dist.kb in upstream.dist.kb.levels) {
  for(conservation in conservation.levels) {
    bp1 = upstream.cons.dist[[ upstream.dist.kb ]]
    upstream.bp =
      apply(bp1[ , colnames(bp1) >= conservation ], 1, sum)
    upstream.bp.at.cutoff[ as.character(upstream.dist.kb),
      as.character(conservation) , ] = upstream.bp[ genes ]
  }
}

save(upstream.bp.at.cutoff,
  file = "git/tf/motif/conservation/upstream.bp.at.cutoff.Rdata")

