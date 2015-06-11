# Filters non-Ce motifs, so as not to be running FIMO on
# all of them.
# (Also omits motifs listed as being Ce motifs.)
# deprecated: now in motifCluster.r

load("git/sort_paper/tf/motif/hughes/motifCluster.Rdata")

# Given a clustering, picks a representative set of motifs.
filter.motifs = function(h) {

  a = cutree(h, h = 0.01)
  a = a[ order(names(a)) ]
  a[ !duplicated(a) ]
}

# Writes out a set of (hopefully) non-redundant motifs.
# Args:
#   org - the organism in question
# Side effects: writes a file ending in "_motif_nr.txt", of
# motifs to search for.
write.nr.motifs = function(org) {
  motif.set = filter.motifs(hughes.motif.cluster[[org]])
  m = setdiff(names(motif.set), hughes.motif.cluster[["Ce"]]$labels)
  write.table(m,
    file=paste0("git/sort_paper/tf/motif/hughes/", org, "_motif_nr.txt"),
    sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# mm.motif.set = filter.motifs(hughes.motif.cluster[["Mm"]])
# write.table(names(mm.motif.set), file="git/sort_paper/tf/motif/hughes/Mm_motif_nr.txt",
#   sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
if (TRUE) {
  write.nr.motifs("Dm")
  write.nr.motifs("Hs")
  write.nr.motifs("Mm")
}

