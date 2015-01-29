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

# first, names from the table of orthologous genes...
r = read.tsv("git/tf/motif.ortholog.3.tsv")
r = r[ , c(2,3) ]
colnames(r) = c("id", "name")
r = r[ !duplicated(r$id) , ]

# also, add in other motif names
r1 = {
  meme.tf.annotate = read.table("git/tf/motif/meme.tf.annotate.tsv",
    sep="\t", header=TRUE, as.is=TRUE)
  r1 = data.frame(id = meme.tf.annotate$id,
    name=paste0(meme.tf.annotate$organism, " ", meme.tf.annotate$gene),
    stringsAsFactors=FALSE)
  r1$name = sub("^ ", "", r1$name)
  r1
}

r = rbind(r, r1)
r = r[ !duplicated(r$id) , ]

motif.name = ifelse(r$name == "", r$id, r$name)
names(motif.name) = r$id

