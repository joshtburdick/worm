# Human-readable names of motifs.
# FIXME: some of these names still look like accession numbers.

source("git/utils.r")

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

motif.ortholog = read.tsv("git/tf/motif.ortholog.3.tsv")

# first, names from the table of orthologous genes...
r = motif.ortholog
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

# For each motif, a list of potential orthologs
orthologs.by.motif = by(motif.ortholog$gene, motif.ortholog$motif.id,
  function(x) {
    x = as.character(x)
    nhrs = grep("nhr", x)
    if (length(x) >= 5 && length(nhrs) >= 5) {
      return(unique(c(grep("nhr", x, invert=TRUE, value=TRUE), paste(length(nhrs), "NHRs"))))
    }

    return(unique(as.character(x))) 
  }
)

