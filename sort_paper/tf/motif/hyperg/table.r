# Writes out results as a table.
# XXX for now, just finding the most significant p-value.

library(reshape2)

source("git/utils.r")

# these are for getting more readable motif names
motif.list = read.table("data/tf/meme/motifList.tsv",
  sep="\t", header=TRUE, as.is=TRUE)
motif.to.name = c(by(motif.list$name,
  motif.list$id,
  function(x) as.character(x)[1]))

motif.ortholog = read.tsv(gzfile("git/tf/motif.ortholog.3.tsv.gz"))
hughes.motif.to.name = c(by(motif.ortholog$motif.name,
  motif.ortholog$motif.id,
  function(x) as.character(x)[1]))

enrich.result.dir = "git/sort_paper/tf/motif/hyperg/allResults/"

# Table including just the most significant p-value.
# XXX possibly not used
enrich.to.table.small = function(enrich) {
  min.p = apply(enrich[,,"p.corr",,,], c(1,2), min)
  melt(min.p)
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




# Given an array of enrichment results, converts it to a table,
# then filters it to only include the most significant.
enrich.to.table = function(enrich) {
  r = dcast(melt(enrich),
    motif + group + upstream.dist + conservation + score ~ stat)
  r = r[ order(r$p) , ]
  r = r[ !duplicated(r[,c("motif","group")]) , ]
  r
}

# Constructs a table of one set of results.
results.table.small = function(name) {

  load(paste0(enrich.result.dir, "5kb/", name, ".Rdata"))
  r1 = enrich.to.table.small(enrich)
  r1$motif.name = motif.to.name[r1$motif]

  load(paste0(enrich.result.dir, "hughes_20141202/", name, ".Rdata"))
  r2 = enrich.to.table.small(enrich)
  r2$motif.name = hughes.motif.to.name[r2$motif]

  r = rbind(r1, r2)
  r = r[ order(r$p) , ]
  r2
#  r1
}

# foo = results.table.small("facs")
load(paste0(enrich.result.dir, "5kb/facs.Rdata"))
enrich = enrich[1:5,1:5,,,,]

