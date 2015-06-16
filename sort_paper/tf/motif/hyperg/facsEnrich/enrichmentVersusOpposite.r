# Computes enrichment of things in e.g. enriched fractions versus depleted,
# comparing, e.g., things enriched versus things depleted.

source("git/tf/motif/enrichment/motifHyperg.r")
source("git/sort_paper/FACS/enrichedInFraction.r")

# Sets of genes enriched, versus depleted.
# Args:
#   r - a data table
#   cutoff - the cutoff to use
# Returns: a list of boolean vectors, one per column of r.
enriched.vs.depleted.1 = function(r, cutoff) {
  a = NULL
  for(s in colnames(r)) {
    a[[ s ]] =
      (r[,s] >= cutoff)[ abs(r[,s]) > cutoff ]
  }
  a
}

gene.sets = {
  plus = enriched.vs.depleted.1(r.sort.only.averaged, 1)
  names(plus) = paste(names(plus), "enriched")
  minus =  enriched.vs.depleted.1(-r.sort.only.averaged, 1)
  names(minus) = paste(names(minus), "depleted")
  c(plus, minus)
}

# XXX hokey way of getting this
orig.motif.list = {
  load("git/sort_paper/tf/allResults/motif/hier.300.clusters.Rdata")
  dimnames(enrich)[[1]]
}

if (FALSE) {
# Looks for motif enrichments for one clustering.
# Args:
#   cl - the clustering
#   name - name to use for the output file
# Side effects: saves motif enrichments to that file.
enrich.test = function(gene.sets, name) {
  output.dir = "git/sort_paper/tf/motif/hyperg/allResults_facs_check/"

  # XXX hack to work around empty gene sets
  num.in.set = sapply(facs.enriched.depleted, sum)
  s = names(num.in.set)[ num.in.set > 0 ]
  gene.sets = gene.sets[ s ]

  system(paste0("mkdir -p ", output.dir, "5kb/"))
  enrich = enrich.test.gene.sets.many.motifs(
    "git/tf/motif/count/upstreamMotifCount/5kb/",
    gene.sets, orig.motif.list)
  save(enrich, file=
    paste0(output.dir, "5kb/", name, ".Rdata"))

  system(paste0("mkdir -p ", output.dir, "hughes/"))
  enrich = enrich.test.gene.sets.many.motifs(
    "git/tf/motif/count/upstreamMotifCount/hughes_20141202/", gene.sets)
  save(enrich, file=
    paste0(output.dir, "hughes/", name, ".Rdata"))

  system(paste0("mkdir -p ", output.dir, "chip/"))
  enrich = enrich.test.gene.sets.many.motifs(
    "git/tf/chip/count/upstreamChipCount/", gene.sets)
  save(enrich, file=
    paste0(output.dir, "chip/", name, ".Rdata"))
}
}
# enrich.test(gene.sets, "facs_enrichmentVsOpposite_2")

