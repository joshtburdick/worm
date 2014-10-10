# Summarizes cases in which a motif was enriched.

require("reshape2")

source("git/utils.r")

# load("git/sort_paper/tf/motifEnrichment/hier.300.clusters.Rdata")

# Given a "slice" of the enrichment results,
# finds the most significant enrichment (by
# corrected p-value), and returns the stats
# about that enrichment.
most.significant.stats.1 = function(a) {

  # find index of most significant enrichment
  i = as.vector(which(a[,,,"p.corr"] == min(a[,,,"p.corr"]),
    arr.ind=TRUE)[1,])

  # converts a cutoff from a number to a numerical level
  cutoff = function(n) as.numeric(dimnames(a)[[n]][ i[n] ])

  c(upstream.dist.kb = cutoff(1),
    conservation = cutoff(2),
    motif.score = cutoff(3),
    a[i[1], i[2], i[3], ])
}

# Given an array of enrichment significances, converts
# to a table containing just the most significant.
# XXX this is a bit hokey, but is at least fast-ish.
most.significant.stats = function(enrich) {

  ms = apply(enrich, c(1,2), most.significant.stats.1)

  r = melt(ms[1,,])[,c(1:2)]
  for(n in dimnames(ms)[[1]]) {
    r = cbind(r, melt(ms[n,,])[,3])
  }
  r[,1] = as.character(r[,1])
  r[,2] = as.character(r[,2])
  colnames(r)[3:ncol(r)] = dimnames(ms)[[1]]
  r
}

# Finds the most significant stats for all the clusters.
most.significant.all.clusterings = function(input.dir, output.dir) {
  system(paste("mkdir -p", output.dir))

  for(f in list.files(input.dir)) {
    write.status(f)
    output.file = sub(".Rdata", ".tsv", f)
    enrich = NULL
    load(paste0(input.dir, "/", f))
    r = most.significant.stats(enrich)
    r = r[ r$p.corr <= 0.05 , ]
    write.tsv(r, paste0(output.dir, "/", output.file))
  }
}

# As above, but only runs on one file, and doesn't do any
# filtering by p-value.
most.significant.all.clusterings.unfiltered = function(input.dir, output.dir, f) {
  system(paste("mkdir -p", output.dir))

  write.status(f)
  output.file = sub(".Rdata", "_unfiltered.tsv.gz", f)
  enrich = NULL
  load(paste0(input.dir, "/", f))
  r = most.significant.stats(enrich)
# not filtering here
  write.tsv(r, gzfile(paste0(output.dir, "/", output.file)))
}

# foo = most.significant.stats(enrich)

if (FALSE) {
most.significant.all.clusterings(
  "git/sort_paper/tf/allResults/motif/",
  "git/sort_paper/tf/summary/motif/")
most.significant.all.clusterings(
  "git/sort_paper/tf/motif/allResults/hughes_motif/",
  "git/sort_paper/tf/motif/summary/hughes_motif/")
}

if (FALSE) {
most.significant.all.clusterings.unfiltered(
  "git/sort_paper/tf/allResults/motif/",
  "git/sort_paper/tf/summary/motif/",
  "hier.300.clusters.Rdata")
}

