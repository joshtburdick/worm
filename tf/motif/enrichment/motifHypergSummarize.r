# Summarizes array of results from a set of motif hypergeometric tests.

source("git/utils.r")
source("git/tf/motif/motifName.r")

library(reshape2)

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
summarize.hyperg.old = function(enrich) {
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


# Given a "slice" of the enrichment results,
# finds the most significant enrichment (by
# corrected p-value), and returns the stats
# about that enrichment.
most.significant.stats.1 = function(a) {

  # find index of most significant enrichment
  i = as.vector(which(a["p",,,] == min(a["p",,,]),
    arr.ind=TRUE)[1,])

  # converts a cutoff from a number to a numerical level
  cutoff = function(n) as.numeric(dimnames(a)[[n+1]][ i[n] ])

  c(upstream.dist.kb = cutoff(1),
    conservation = cutoff(2),
    motif.score = cutoff(3),
    a[ , i[1], i[2], i[3] ])
}

# Given an array of enrichment results, gets information about
# the most significant result.
enrich.to.table = function(enrich) {
  enrich1 = apply(enrich, c(1,2), most.significant.stats.1)
  r = dcast(melt(enrich1), motif + group ~ Var1)
  r$motif = as.character(r$motif)
  r$group = as.character(r$group)
#  r$motif.name = motif.name[ r$motif ]
#  r$motif.name[ is.na(r$motif.name) ] = ""
#  r = r[ , c(1,12,2:11) ]
  r
}

# Same as above, but combines several arrays' worth of results,
# correcting p-values by fdr together.
# XXX this is slow, and uses a lot of memory (particularly
# "enrich.to.table()").
enrich.to.table.many = function(enrich.list) {

  p = sapply(enrich.list, function(e) as.vector(e[,,"p",,,]))

  # XXX compute indices
  a = cumsum(c(1, sapply(p, length)))
  b = cumsum(sapply(p, length))

  # adjust all p values together
  p = as.numeric(c(p, recursive=TRUE))
cat("concatted p")
  p = p.adjust(p, method="fdr")
cat("adjusted p")

  r = NULL
  for(i in 1:length(enrich.list)) {
cat("converting", names(enrich.list)[[i]])
    e1 = enrich.list[[i]]
cat("created array")
    e1[,,"p.corr",,,] = p[ a[i] : b[i] ]
cat("wrote to array")
    t1 = enrich.to.table(e1)
cat("converted to table")
    r = rbind(r, t1)
cat("did rbind")
  }

  r
}

# Same as above, but combines several arrays' worth of results,
# correcting p-values by fdr together.
# XXX this is slow, and uses a lot of memory
enrich.to.table.many.1 = function(enrich.list) {

  p = sapply(enrich.list, function(e) as.vector(e[,,"p",,,]))
  p1 = c(p, recursive=TRUE)
  p2 = p.adjust(p1, method="fdr")

  # XXX compute indices
  a = cumsum(c(1, sapply(p, length)))
  b = cumsum(sapply(p, length))
  r = NULL
  for(i in 1:length(enrich.list)) {
    enrich.list[[i]][,,"p.corr",,,] = p2[ a[i] : b[i] ]
    r = rbind(r, enrich.to.table(enrich.list[[i]]))
  }

  r
}
