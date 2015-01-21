# Combines the original motifs and Hughes motifs into one
# file.
# (For now, only doing this for the 300-cluster version.)

source("git/utils.r")
source("git/sort_paper/tf/most.significant.r")

motif.ortholog = read.tsv("git/tf/motif.ortholog.3.tsv")

system("mkdir -p git/sort_paper/tf/motif/combined/")

# Adjusts p-values. Any p-values less than 1e-300 first get
# rewritten to 1e-300.
p.adjust.1 = function(p) {
  p[ is.na(p) ] = 1
  p[ p <= 1e-300 ] = 1e-300
  p.adjust(p, method="fdr")
}

# Combines one pair of clusterings.
# combine = function(clustering.name) {
clustering.name = "hier.300.clusters"

  # the original clustering
  load(paste0("git/sort_paper/tf/allResults/motif/",
    clustering.name, ".Rdata"))
  motif.enrich.1 = enrich

  # results using the Hughes motifs
  load(paste0("git/sort_paper/tf/motif/allResults/hughes_motif/",
    clustering.name, ".Rdata"))
  hughes.enrich = enrich

  # limit Hughes motif to those for which we have some annotation
  m = intersect(motif.ortholog[,"motif.id"], dimnames(hughes.enrich)[[1]])
  hughes.enrich = hughes.enrich[m,,,,,]

  # re-correct p-values
  motif.enrich.1[,,,,,"p.corr"] = p.adjust.1(motif.enrich.1[,,,,,"p"])
  hughes.enrich[,,,,,"p.corr"] = p.adjust.1(hughes.enrich[,,,,,"p"])

  # FIXME: write out the "appended" intermediate file


  # compute most significant
  motif.enrich.ms = most.significant.stats(motif.enrich.1)
  hughes.enrich.ms = most.significant.stats(hughes.enrich)

  write.tsv(rbind(motif.enrich.ms, hughes.enrich.ms), gzfile(
    paste0("git/sort_paper/tf/motif/combined/", clustering.name, ".tsv")))

# }


# combine("hier.300.clusters")

