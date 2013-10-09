# Constructs a network based on enriched upstream motifs.

source("git/utils.r")



# Reads in the various enrichments, and writes out a file
# which Cytoscape can read.
# Args:
#   path - where the various enrichment files are
# Side effects: writes a Cytoscape-readable file in that directory
write.enrichment.cytoscape.file = function(path) {

  r1 = read.tsv(paste(path, "chipEnrichment_5kb.tsv", sep="/"))
  r1$type = "ChIP"
  r2 = read.tsv(paste(path, "deNovoMotifEnrichment_5kb.tsv", sep="/"))
  r2$type = "deNovoMotif"
  r3 = read.tsv(paste(path, "knownMotifEnrichment_5kb.tsv", sep="/"))
  r3$type = "knownMotif"

  r = rbind(r3)       # XXX temporary
  r = r[ r$p.bh <= 0.05, ]
  r1 = r[ , c("name", "type", "cluster") ]
  write.tsv(r1, paste(path, "cytoscape.tsv", sep="/"))
}


write.enrichment.cytoscape.file(
  "git/cluster/hierarchical/hier.ts.300.clusters/")

