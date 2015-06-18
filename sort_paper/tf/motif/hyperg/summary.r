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
# This also filters motifs to be non-redundant (which shouldn't
# be neccessary, admittedly.)
hughes.enrichments = function(enrich.dir, clustering) {
  enrich.list = list()
  m = c()   # motifs seen so far
#  for(org in c("Mm", "Dm")) {
  for(org in c("Ce", "Dm", "Hs", "Mm")) {
    motif.set = paste0(org, "_1.02")
    write.status(paste(clustering, motif.set))
    enrich = NULL
    load(paste0(enrich.dir, "/", motif.set, "/", clustering, ".Rdata"))
# reduce the non-Ce motifs, based on motif clustering,
# which arguably should have been done earlier.
# also omit motifs we've seen before.
    if (org %in% names(nr.motifs)) {
      m1 = intersect(dimnames(enrich)[[1]], nr.motifs[[org]])
      m1 = setdiff(m1, m)
      enrich = enrich[m1,,,,,]
      m = union(m, m1)
    }
    enrich.list[[motif.set]] = enrich
  }
  enrich.list
}


if (TRUE) {
  enrich.list = hughes.enrichments(enrich.dir, "hier.300.clusters")
  r = enrich.to.table.many.1(enrich.list)
  write.tsv(r, file=gzfile(paste0(output.dir, "hughes/hier.300.clusters.tsv.gz")))
}

if (FALSE) {
  system(paste0("mkdir -p ", output.dir, "/hughes/"))
  for(cl in c("facs_vs_opposite_1")) {
#  for(cl in c("spencer", "spencerEmbryonic", "hy")) {
    enrich.list = hughes.enrichments(enrich.dir, cl)
    r = enrich.to.table.many.1(enrich.list)
    write.tsv(r,
      file=gzfile(paste0(output.dir, "/hughes/", cl, ".tsv.gz")))
  }
}



