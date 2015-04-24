# Summarizes array of results from a set of motif hypergeometric tests.

source("git/utils.r")
source("git/tf/motif/motifName.r")

# this is deprecated
if (FALSE) {
hughes.motif.info = read.table(
  "data/tf/hughes/Caenorhabditis_elegans_2014_12_02_11_30_am/TF_Information_all_motifs.txt.gz",
  sep="\t", header=TRUE, as.is=TRUE)

motif.id.2.name = c(by(hughes.motif.info$DBID.1,
  hughes.motif.info$Motif_ID,
  function(x) as.character(x)[1]))
}

# Summarizes most significant results from an enrichment test.
# Args:
#   enrich - the results from, e.g., enrich.test.gene.sets.many.motifs()
summarize.hyperg = function(enrich) {
  # ??? not clear if "drop=FALSE" is always correct here
  min.p =
    apply(enrich[,,"p.corr",,,,drop=FALSE], c(1,2), min)

  r = data.frame(
    motif = rep(rownames(min.p), times=ncol(min.p)),
    group = rep(colnames(min.p), each=nrow(min.p)),
    p.corr = as.vector(min.p),
    stringsAsFactors = FALSE)
  r$motif.name = motif.name[ r$motif ]
  r$motif.name[ is.na(r$motif.name) ] = ""
  r = r[ order(r$p.corr),  ]   # was c(2,1,4,3)
}


