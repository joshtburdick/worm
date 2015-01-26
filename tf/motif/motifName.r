# Human-readable names of motifs.

# FIXME: some of these names still look like accession numbers.

# Human-readable motif names (old).
if (FALSE) {
motif.name = {
  r = read.table("git/tf/motif/motifList.tsv",
    header=TRUE, fill=TRUE, as.is=TRUE)
  r = r[,c(2:3)]
  r = r[ !duplicated(r$id) , ]
  motif.name = ifelse(r$name == "", r$id, r$name)
  names(motif.name) = r$id
  motif.name
}
}

r = read.tsv("git/tf/motif.ortholog.3.tsv")
r = r[ , c(2,3) ]
colnames(r) = c("id", "name")
r = r[ !duplicated(r$id) , ]

motif.name = ifelse(r$name == "", r$id, r$name)
names(motif.name) = r$id

