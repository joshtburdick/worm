# Summarizes cases in which a motif was enriched.

require("reshape2")

source("git/utils.r")
source("git/util/arrayUtils.r")


source("git/sort_paper/tf/most.significant.r")

# load("git/sort_paper/tf/motifEnrichment/hier.300.clusters.Rdata")


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

# Similar, but includes "known" set and Hughes motifs.
# Doesn't filter by significance.
# Args:
#   output.dir - where to write the results
#   f - the clustering to do this for
# Side effects - writes a .tsv.gz file in output.dir
most.significant.merged = function(output.dir, f) {
  system(paste("mkdir -p", output.dir))
  write.status(f)
  output.file = sub(".Rdata", "_merged.tsv.gz", f)
  enrich = NULL
  load(paste0("git/sort_paper/tf/allResults/motif/", f))
  enrich.1 = enrich
  enrich = NULL
  load(paste0("git/sort_paper/tf/motif/allResults/hughes_motif/", f))
  enrich.2 = enrich
  enrich = array.bind(list(enrich.1, enrich.2), 1)
  enrich[,,,,,"p.corr"] = array(p.adjust(as.vector(enrich[,,,,,"p"]),
    method="fdr"),
    dim=dim(enrich[,,,,,"p.corr"]),
    dimnames=dimnames(enrich[,,,,,"p.corr"]))

  r = most.significant.stats(enrich)
  colnames(r)[[1]] = "motif"
  colnames(r)[[2]] = "group"

# not filtering by significance here
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

if (TRUE) {
most.significant.merged(
  "git/sort_paper/tf/summary/motif/",
  "hier.300.clusters.Rdata")
}
