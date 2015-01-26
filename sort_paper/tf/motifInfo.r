# Utilities for getting information about motifs.

source("git/utils.r")

# information about orthologs
#motif.ortholog = read.table("git/tf/motif.ortholog.3.tsv",
#  as.is=TRUE, header=TRUE)
#
motif.ortholog = read.tsv("git/tf/motif.ortholog.3.tsv")

motif.filter = read.tsv("git/tf/motif/motifFilter.tsv")

# add representative motif name, from the motif clustering
motif.ortholog$canonical.motif =
  motif.filter[ motif.ortholog$motif, "canonical.name" ]
motif.ortholog = motif.ortholog[ !is.na(motif.ortholog$canonical.motif) , ]
# for each motif, list of potential orthologs
motif.info = by(motif.ortholog$gene, motif.ortholog$canonical.motif,
  function(g) {
    g = unique(as.character(g))

    # summarize NHRs
    nhrs = grep("nhr", g)
    if (length(g) >= 5 && length(nhrs) >= 5) {
      return(unique(c(grep("nhr", g, invert=TRUE, value=TRUE), paste(length(nhrs), "NHRs"))))
    }

    # if there are many genes, summarize this
    if (length(g) > 5) {
      g = c(g[1:4], paste("and", length(g)-4, "others"))
    }

    return(paste(g, collapse=" "))
  }
)
motif.info = c(motif.info, recursive=TRUE)

# The gene for each motif.
meme.tf.annotate = read.table("git/tf/motif/meme.tf.annotate.tsv",
  sep="\t", header=TRUE, as.is=TRUE)
motif.gene = by(meme.tf.annotate$gene, meme.tf.annotate$id,
  function(x) as.character(x[1]))

