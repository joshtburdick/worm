# Computes counts of ChIP peaks upstream of genes, at
# different conservation and distance cutoffs.
# XXX I'm sometimes referring to "ChIP experiments" as
# "motifs", which is somewhat abuse of notation.

source("git/utils.r")

# list of gene names to include
gene.names = read.table(
  "git/tf/motif/motifCount/regions/upstream_5kb_WS220.bed",
  sep="\t", as.is=TRUE)[,4]

motif.gene.dir = "git/cluster/chip/distAndConservation new/"

# Gets motif occurrences relative to genes from a file
# (the output of "distAndConservation.py").
get.motifs.from.file = function(motif.gene.dir) function(motif) {
  r = read.table(gzfile(
    paste(motif.gene.dir, motif, "_upstreamChipCons.tsv.gz", sep="")),
    as.is=TRUE)
  colnames(r) =
    c("region.chr", "region.a", "region.b",
    "gene", "score", "strand",
    "motif", "motif.a", "motif.b", "motif.id", "motif.score",
    "motif.strand", "motif.cons")
  r$upstream.dist = ifelse(r$strand=="+",
    r$motif.b - r$region.b, r$region.a - r$motif.a)
  r
}

# Counts motifs at particular cutoffs.
# Args:
#   motifs - the motifs to include
#   gene.names - names of genes to include
#   cutoffs - named list of vectors, including fields
#     "upstream.dist", "conservation", and "score",
#      giving cutoffs to use for those
#   get.motifs.f - function which, given a motif name,
#     gets the counts for that motif
#   output.dir - where to write output
# Side effects: writes out motif counts in the output directory
count.motifs.at.cutoffs =
    function(motifs, gene.names, cutoffs, get.motifs.f, output.dir) {

  system(paste("mkdir -p", output.dir))

  # loop through the motifs and cutoffs
  for (motif in motifs) {

    # array of results
    motif.count = array.from.dimnames(list(gene = gene.names,
      upstream.dist = cutoffs$upstream.dist,
      conservation = cutoffs$conservation,
      score = cutoffs$score))

    m = get.motifs.f(motif)
    for (upstream.dist.kb in cutoffs$upstream.dist)
      for (conservation in cutoffs$conservation)
        for (motif.score in cutoffs$score) {
          write.status(paste(motif, upstream.dist.kb,
            conservation, motif.score))
          upstream.dist.bp = -1000 * upstream.dist.kb
          m1 = m[ m$upstream.dist >= upstream.dist.bp &
            m$motif.cons >= conservation , "gene" ]
          mc1 = c(table(m1))
          g = intersect(names(mc1), gene.names)
          motif.count[ , as.character(upstream.dist.kb),
            as.character(conservation),
            as.character(motif.score) ] = 0
          motif.count[ g, as.character(upstream.dist.kb),
            as.character(conservation),
            as.character(motif.score) ] = mc1[g]
          gc()    # this seems faster
        }

      save(motif.count, file=paste0(output.dir, "/", motif, ".Rdata"),
        compress="bzip2")
    }
}

chip.experiments =
  sub("_upstreamChipCons.tsv.gz", "", list.files(motif.gene.dir))

if (TRUE) {
r = count.motifs.at.cutoffs(chip.experiments,
  gene.names,
  list(upstream.dist = c(1, 2, 3),
    conservation = c(0, 0.5, 0.7, 0.9),
    score = c(0)),
  get.motifs.from.file(motif.gene.dir),
  "git/tf/chip/count/upstreamChipCount/")
}

