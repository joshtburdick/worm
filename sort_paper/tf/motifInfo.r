# Utilities for getting information about motifs.

source("git/utils.r")

# previous version of motif-ortholog stuff
if (FALSE) {
  # information about orthologs
  motif.ortholog = read.tsv("git/tf/motif.ortholog.3.tsv")
  motif.filter = read.tsv("git/tf/motif/motifFilter.tsv")

  # add representative motif name, from the motif clustering
  motif.ortholog$canonical.motif =
    motif.filter[ motif.ortholog$motif.id, "canonical.name" ]
  motif.ortholog = motif.ortholog[ !is.na(motif.ortholog$canonical.motif) , ]
}

source("git/sort_paper/tf/motif/hughes/motifInfo.r")

# For each motif, a concise list of potential orthologs.
ortholog.by.motif.small = by(motif.ortholog$gene, motif.ortholog$motif.id,
  function(g) {
    g = unique(as.character(g))

    # slight reordering
    g = g[ order( !(g %in% c("pha-4"))) ]

    # summarize NHRs
#    nhrs = grep("nhr", g)
#    if (length(g) >= 5 && length(nhrs) >= 5) {
#      return(unique(c(grep("nhr", g, invert=TRUE, value=TRUE), paste(length(nhrs), "NHRs"))))
#    }

    # if there are many genes, summarize this
    if (length(g) > 5) {
      g = c(g[1:4], paste("and", length(g)-4, "others"))
    }

    return(paste(g, collapse=" "))
  }
)
ortholog.by.motif.small = c(ortholog.by.motif.small, recursive=TRUE)

# The gene for each motif.
motif.gene = c(by(motif.info.1$related.gene, motif.info.1$motif.id,
  function(x) as.character(x[1])))

