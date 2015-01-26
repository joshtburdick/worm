# Computes motif counts.
# Currently "countMotifs.pl" does most of the work.

source("git/utils.r")
source("git/data/name_convert.r")

load("git/tf/motif/clusterUsingMotIV.Rdata")

# Processes the output of the motif search slightly.
process.motif.file = function(f, region.sizes) {
#  m = t( read.table(f, sep="\t", as.is=TRUE, check.names=FALSE) )
  m = t(read.tsv(f))

  # rename genes to match up with the clustered expression data
  # (ideally the gene locations should have those names anyway)
  m = rename.gene.names(m)

  # normalize this
  m = m / region.sizes

  # ??? only keep motifs which we clustered?
  m[ , intersect(colnames(m), motif.clusters$label) ]
#  m
}

# Gets sizes of exons in a BED file.
# Args:
#   f - a BED file
# Returns: sizes of genes
bed.gene.sizes = function(f) {
  r = read.table(f, sep="\t", as.is=TRUE)

  # if there are only six columns, use the bounds of
  # the entire region
  if (ncol(r) == 6) {
    size = data.frame(s = r[,3] - r[,2], foo=1, row.names=r[,4])
    size = rename.gene.names(size)
    s = size$s
    names(s) = rownames(size)
    return(s)
  }

  # otherwise, total up the exon sizes
  if (ncol(r) == 12) {
    size = data.frame(s =
      sapply(strsplit(r[,11], ","), function(x) sum(as.numeric(x))), foo=1,
      row.names = r[,4])
    size = rename.gene.names(size)
    s = size$s
    names(s) = rownames(size)
    return(s)
  }

  return(NULL)
}

# get sizes of upstream regions
upstream.region.size = bed.gene.sizes(
  "git/tf/motif/motifCount/upstream_liftOver_WS220.bed")
upstream.region.size.0.5.cons = bed.gene.sizes(
  "git/tf/motif/motifCount/upstream_liftOver_WS220_0.5cons.bed")

if (FALSE) {
known.motifs = process.motif.file(gzfile(
  "git/tf/motif/motifCount/knownMotif_5kbUp.tsv.gz"),
  upstream.region.size)
de.novo.motifs = process.motif.file(gzfile(
  "git/tf/motif/motifCount/deNovoMotif_5kbUp.tsv.gz"),
  upstream.region.size)

known.motifs.0.5.cons = process.motif.file(gzfile(
  "git/tf/motif/motifCount/knownMotif_5kbUp_0.5cons.tsv.gz"),
  upstream.region.size.0.5.cons)
de.novo.motifs.0.5.cons = process.motif.file(gzfile(
  "git/tf/motif/motifCount/deNovoMotif_5kbUp_0.5cons.tsv.gz"),
  upstream.region.size.0.5.cons)

save(known.motifs, de.novo.motifs,
  known.motifs.0.5.cons, de.novo.motifs.0.5.cons,
  upstream.region.size, upstream.region.size.0.5.cons,
  file="git/tf/motif/motifCount/motif.counts.Rdata")
}


if (FALSE) {
shuffled.known.motifs = process.motif.file(gzfile(
  "git/tf/motif/motifCount/shuffledKnownMotif_5kbUp.tsv.gz"),
  upstream.region.size)
shuffled.known.motifs.0.5.cons = process.motif.file(gzfile(
  "git/tf/motif/motifCount/shuffledKnownMotif_5kbUp_0.5cons.tsv.gz"),
  upstream.region.size.0.5.cons)

save(shuffled.known.motifs, shuffled.known.motifs.0.5.cons,
  upstream.region.size, upstream.region.size.0.5.cons,
  file="git/tf/motif/motifCount/shuffled.motif.counts.Rdata")
}



