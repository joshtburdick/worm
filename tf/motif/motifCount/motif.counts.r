# Computes motif counts.
# Currently "countMotifs.pl" does most of the work.

load("git/tf/motif/clusterUsingMotIV.Rdata")

known.motifs = {

  # the known motifs, 5kb upstream (without conservation filtering)
  m = read.table(gzfile("git/tf/motif/motifCount/knownMotif_5kbUp.tsv.gz"),
    sep="\t", header=TRUE, row.names=1, check.names=FALSE, fill=TRUE)

  # omit cases in which motif search didn't work for all genes
  # XXX for whatever reason, some motifs such as MX000045 are truncated.
  # I have no idea why this is, but as it's part of PRODORIC
  # (a prokaryotic database) anyway, for now I'm just ignoring the
  # problem (it will be omitted anyway.) (The issue may be
  # because that motif occurs a lot.)
  m = m[ apply(is.na(m), 1, sum) == 0 , ]

  # only include motifs which were clustered
  m = m[ intersect(rownames(m), motif.clusters$label) , ]

  t(m)
}

de.novo.motifs = {
  m = read.table(
    gzfile("git/tf/motif/motifCount/deNovoMotif_5kbUp.tsv.gz"),
    sep="\t", header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE)

  # see above
  m = m[ apply(is.na(m), 1, sum) == 0 , ]

  t(m)
}

save(known.motifs, de.novo.motifs, file="git/tf/motif/motifCount/motif.counts.Rdata")

