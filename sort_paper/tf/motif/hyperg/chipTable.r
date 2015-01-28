# Writes out ChIP enrichment results as a table.
# XXX this is fairly similar to the other table, but somewhat
# different, as the ChIP data is a somewhat different shape

library(reshape2)

source("git/utils.r")

enrich.result.dir = "git/sort_paper/tf/motif/hyperg/allResults/"

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
    a[ , i[1], i[2], 1 ])
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

most.significant.results.chip = function(name) {

  load(paste0(enrich.result.dir, "chip/", name, ".Rdata"))
  r = enrich.to.table(enrich)

  colnames(r)[1] = "chip experiment"
  colnames(r)[6] = "genes with peak in cluster"
  colnames(r)[7] = "genes in cluster"
  colnames(r)[8] = "genes with peak"
  r = r[ order(r$p) , c(1:4,6:10) ]
  r
}

r = most.significant.results.chip("hier.300")
system("mkdir -p git/sort_paper/tf/motif/hyperg/chipTable")
write.tsv(r, gzfile("git/sort_paper/tf/motif/hyperg/chipTable/hier.300.tsv.gz"))


