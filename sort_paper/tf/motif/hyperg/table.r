# Writes out results as a table.
# XXX for now, just finding the most significant p-value.

library(reshape2)

source("git/utils.r")

# for getting more readable motif names
motif.list = read.table("data/tf/meme/motifList.tsv",
  sep="\t", header=TRUE, as.is=TRUE)
motif.to.name.old = c(by(motif.list$name,
  motif.list$id,
  function(x) as.character(x)[1]))

motif.ortholog = read.tsv("git/tf/motif.ortholog.3.tsv")
motif.to.name = c(by(motif.ortholog$motif.name,
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



# Given an array of enrichment results, gets information about
# the most significant result.
enrich.to.table = function(enrich) {
  enrich1 = apply(enrich, c(1,2), most.significant.stats.1)
  r = dcast(melt(enrich1), motif + group ~ Var1)
  r$motif = as.character(r$motif)
  r$group = as.character(r$group)
  r
}

# Constructs a table of one set of results.
most.significant.results = function(name) {

  load(paste0(enrich.result.dir, "5kb/", name, ".Rdata"))
  r1 = enrich.to.table(enrich)
  load(paste0(enrich.result.dir, "hughes/", name, ".Rdata"))
  r2 = enrich.to.table(enrich)
  r = rbind(r1, r2)

  # hopefully improve motif names
  r$motif.name = r$motif
  i = r$motif %in% names(motif.to.name)
  r[ i , "motif.name" ] = motif.to.name[ r[ i , "motif" ] ]

  colnames(r)[6] = "genes with motif in cluster"
  colnames(r)[7] = "genes in cluster"
  colnames(r)[8] = "genes with motif"
  r = r[ order(r$p) , c(1,11,2:10) ]
  r
}

# r = most.significant.results("hier.300")
# system("mkdir -p git/sort_paper/tf/motif/hyperg/table/")
# write.tsv(r, gzfile("git/sort_paper/tf/motif/hyperg/table/hier.300.tsv.gz"))

r = most.significant.results("facs")
write.tsv(r, gzfile("git/sort_paper/tf/motif/hyperg/table/facs.tsv.gz"))


