# Filters non-Ce motifs, so as not to be running FIMO on
# all of them.

load("git/sort_paper/tf/motif/hughes/motifCluster.Rdata")

# Given a clustering, picks a representative set of motifs.
filter.motifs = function(h) {

  a = cutree(h, h = 0.01)
  a = a[ order(names(a)) ]
  a[ !duplicated(a) ]
}


hs.motif.set = filter.motifs(hughes.motif.cluster[["Hs"]])
write.table(names(hs.motif.set), file="git/sort_paper/tf/motif/hughes/Hs_motif_nr.txt",
  sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

mm.motif.set = filter.motifs(hughes.motif.cluster[["Mm"]])
write.table(names(mm.motif.set), file="git/sort_paper/tf/motif/hughes/Mm_motif_nr.txt",
  sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

