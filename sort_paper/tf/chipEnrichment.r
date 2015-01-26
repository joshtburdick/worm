# Looks for ChIP enrichment upstream of various sets
# of genes.

source("git/sort_paper/tf/upstreamEnrichment.r")
source("git/sort_paper/FACS/enrichedInFraction.r")

source("git/sort_paper/enrichment/groupToList.r")

clustering.dir = "git/cluster/hierarchical/"

chip.dir = "git/cluster/chip/distAndConservation/"
# chip.dir = "git/cluster/chip/distAndConservation_HOT/"

chip.files = sort(sub("_upstreamChipCons.tsv.gz", "",
  list.files(chip.dir)))
# XXX removing some duplicated filenames; probably I
# should just remove the underlying files
chip.files = setdiff(chip.files,
  c("AMA-1_L4-Young-Adult-larvae_rep1",
    "AMA-1_Young-adult_rep1",
    "DAF-16_L4-Young-Adult-larvae_rep1",
    "DPL-1_Young-adult_rep1",
    "DPL-1_Young-adult_rep2",
    "EFL-1_Young-adult_rep1",
    "EFL-1_Young-adult_rep2",
    "PHA-4_Young-adult_rep1",
    "PHA-4_Young-adult_rep2",
    "W03F9.2_L4-Young-Adult-larvae_rep1",
    "W03F9.2_L4-Young-Adult-larvae_rep2"))

chip.file.sets = by(chip.files, sub("_rep\\d", "", chip.files),
  function(x) c(as.character(x)))
chip.names = names(chip.file.sets)

# XXX
# chip.names = setdiff(chip.names,
#   c("DPL-1_Young-Adult", "EFL-1_Young-Adult", "PHA-4_Young-Adult",
#     "W03F9.2_L4-Young-Adult-stage-larvae"))


# Gets the ChIP peaks upstream of all genes, for one file.
get.chip.locs.one.file = function(m) {

  r = read.table(paste(chip.dir, m, "_upstreamChipCons.tsv.gz", sep=""),
    as.is=TRUE)
  colnames(r) =
    c("region.chr", "region.a", "region.b",
    "gene", "score", "strand",
    "motif", "motif.a", "motif.b", "motif.id", "motif.score",
    "motif.strand", "motif.cons")
  r$upstream.dist = ifelse(r$strand=="+",
    r$motif.b - r$region.b, r$region.a - r$motif.a)
  r
}

# Combines all motifs from all the files.
get.chip.locs = function(m) {
  r = NULL
  for (f in chip.file.sets[[m]]) {
    r = rbind(r, get.chip.locs.one.file(f))
  }

  r
}

# output.dir = "git/sort_paper/tf/allResults/chip/"
output.dir = "git/sort_paper/tf/allResults/chip_incl_HOT/"

system(paste("mkdir -p", output.dir))

if (TRUE) {
e = compute.enrichment.diff.cutoffs(chip.names, get.chip.locs,
  facs.enriched.depleted, c(1,2,3), c(0, 0.5, 0.7, 0.9), c(0))
e[,,,,,"p"][ is.na(e[,,,,,"p"]) ] = 1
e[,,,,,"p.corr"] = array(p.adjust(as.vector(e[,,,,,"p"]), method="fdr"),
    dim=dim(e[,,,,,"p.corr"]), dimnames=dimnames(e[,,,,,"p.corr"]))
enrich = e
save(enrich, file="git/sort_paper/tf/allResults/chip/facs.Rdata")

}

if (TRUE) {

for (f in list.files(clustering.dir)) {
  gc()
  cat(paste(f, date()), "\n")

  # clustering to use
#  clustering1 = read.tsv(paste0(clustering.dir, "/", f, "/clusters.tsv"))
#  clustering = clustering1[,2]
#  names(clustering) = rownames(clustering1)
  clustering = cluster.to.gene.list(paste0(clustering.dir, "/", f, "/clusters.tsv"))

  enrich = compute.enrichment.diff.cutoffs(chip.names, get.chip.locs,
    clustering, c(1,2,3), c(0, 0.5, 0.7, 0.9), c(0))
  enrich[,,,,,"p"][ is.na(enrich[,,,,,"p"]) ] = 1

  enrich[,,,,,"p.corr"] = array(p.adjust(as.vector(enrich[,,,,,"p"]), method="fdr"),
    dim=dim(enrich[,,,,,"p.corr"]), dimnames=dimnames(enrich[,,,,,"p.corr"]))

  save(enrich, file=paste0(output.dir, "/", f, ".Rdata"))
}

}

