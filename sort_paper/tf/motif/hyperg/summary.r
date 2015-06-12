# Summarizes motifs enriched upstream of clusters.

source("git/util/arrayUtils.r")
source("git/tf/motif/enrichment/motifHypergSummarize.r")
source("git/sort_paper/tf/motif/hughes/motifInfo.r")

load("git/sort_paper/tf/motif/hughes/motifCluster.Rdata")

enrich.dir = "git/sort_paper/tf/motif/hyperg/allResults2/"
output.dir = "git/sort_paper/tf/motif/hyperg/summary/"

system(paste0("mkdir -p ", output.dir))

clustering.name = "hier.300.clusters"

# Gets all the Hughes motif enrichment.

hughes.enrichments = function(enrich.dir, clustering) {
  enrich.list = list()
#  for(motif.set in paste0(c("Mm", "Dm"), "_1.02")) {
  for(org in c("Dm", "Hs", "Mm", "Ce")) {
    motif.set = paste0(org, "_1.02")
    write.status(paste(clustering, motif.set))
    enrich = NULL
    load(paste0(enrich.dir, "/", motif.set, "/", clustering, ".Rdata"))
# reduce the non-Ce motifs, based on motif clustering,
# which arguably should have been done earlier.
    if (org %in% names(nr.motifs)) {
      m1 = intersect(dimnames(enrich)[[1]], nr.motifs[[org]])
      enrich = enrich[m1,,,,,]
    }
    enrich.list[[motif.set]] = enrich
  }
  enrich.list
}

enrich.list = hughes.enrichments(enrich.dir, "hier.300.clusters")

r = enrich.to.table.many.1(enrich.list)
write.tsv(r, file=gzfile(paste0(output.dir, "hier.300.clusters.tsv.gz")))


