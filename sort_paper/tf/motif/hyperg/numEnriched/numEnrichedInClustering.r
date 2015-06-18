# Number of motifs enriched upstream of different clusterings.

source("git/utils.r")

enrich.dir =
  "git/sort_paper/tf/motif/hyperg/numEnriched/clusteringComparison/"

r = NULL
for(f in setdiff(list.files(enrich.dir), "Ce_1.02")) {   # XXX workaround
  cl.name = sub(".tsv.gz", "", f)
write.status(f)
  a = read.tsv(paste0(enrich.dir, f))

  for(p.cutoff in c(0.05, 1e-2, 1e-3, 1e-5, 1e-10)) {
    a1 = a[ a$p.corr <= p.cutoff , ]
    r1 = c(clustering.name = cl.name,
      p.cutoff = p.cutoff,
      num.pairs = nrow(a1),
      num.motifs = length(unique(a1$motif)),
      num.clusters = length(unique(a1$group)))
    r = rbind(r, r1)
  }

}

write.tsv(r, "git/sort_paper/tf/motif/hyperg/numEnriched/numEnrichedInClustering.tsv")






