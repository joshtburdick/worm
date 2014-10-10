# Human-readable names of motifs.

# Human-readable motif names.
motif.name = {
  r = read.table("git/tf/motif/motifList.tsv",
    header=TRUE, fill=TRUE, as.is=TRUE)
  r = r[,c(2:3)]
  r = r[ !duplicated(r$id) , ]
  motif.name = ifelse(r$name == "", r$id, r$name)
  names(motif.name) = r$id
  motif.name
}

