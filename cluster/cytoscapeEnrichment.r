# Constructs a network based on things (motifs and ChIP signal)
# upstream of genes in clusters, filtered by correlation of a TF
# with clusters.

source("git/utils.r")

r = read.tsv("git/cluster/readRatios.tsv")
ortho = read.table("git/tf/motif.ortholog.2.tsv", header=TRUE, as.is=TRUE)

# for now, not including InParanoid orthology
ortho = ortho[ ortho$method != "InParanoid" , ]

# XXX some of these motifs aren't included (I'm not sure why),
# so omit those cases
ortho = ortho[ ortho$motif %in% colnames(known.motifs) , ]

# just keep cases in which there's a motif, and only keep
# the motif from the "best" homolog (by "gene.score"; arguably
# these should be normalized in some way)
ortho = ortho[
  order(ortho$gene,
    ortho$species,
    -ortho$related.gene.score) , ]

# ortho = ortho[ !duplicated(ortho[ , c("gene", "species"),]) , ]
# ??? just keeping unique gene names
ortho = ortho[ !duplicated(ortho[ , c("gene"),]) , ]

# the ratio data
r = read.tsv("git/cluster/readRatios.tsv")

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


if (FALSE) {
write.enrichment.cytoscape.file(
  "git/cluster/hierarchical/hier.ts.300.clusters/")
}


