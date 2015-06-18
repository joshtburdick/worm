# Does the hypergeometric test for some clusterings
# (largely to compare how many things were enriched.)

source("git/utils.r")
source("git/tf/motif/enrichment/motifHyperg.r")
source("git/tf/motif/enrichment/motifHypergSummarize.r")

output.dir =
  "git/sort_paper/tf/motif/hyperg/numEnriched/clusteringComparison/"
system(paste("mkdir -p ", output.dir))

# Does enrichment tests for one clustering, and one set of motifs.
hyperg.test.1 = function(cl, motif.set.name, clustering.name) {
  motif.count.base =
    "/home/jburdick/gcb/git/sort_paper/tf/motif/upstreamCount/"
  motif.count.dir = paste0(motif.count.base, "/", motif.set.name, "/")

  system(paste0("mkdir -p ", output.dir, "/"))
  enrich = enrich.test.many.motifs(motif.count.dir, cl)

  r = enrich.to.table(enrich)
  write.tsv(r, gzfile(paste0(output.dir, clustering.name, ".tsv.gz")))
}

# Run this for all clusterings (just for the Ce motifs), to
# see which has the most enrichments.
if (TRUE) {
  for(clustering.name in list.files("git/cluster/hierarchical")) {
    cat(clustering.name, "\n")
    cl = as.matrix(read.tsv(
      paste0("git/cluster/hierarchical/", clustering.name, "/clusters.tsv")))[,2]
    cl = gsub(" ", "", cl)
    hyperg.test.1(cl, "Ce_1.02", clustering.name)
  }
}


