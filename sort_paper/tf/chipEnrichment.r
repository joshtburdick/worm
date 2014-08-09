# Looks for ChIP enrichment upstream of various sets
# of genes.

source("git/sort_paper/tf/upstreamEnrichment.r")
source("git/sort_paper/FACS/enrichedInFraction.r")

source("git/sort_paper/enrichment/groupToList.r")

clustering.dir = "git/cluster/hierarchical/"

chip.dir = "git/cluster/chip/distAndConservation/"
chip.names = sub("_upstreamChipCons.tsv.gz", "",
  list.files(chip.dir))
# XXX
chip.names = setdiff(chip.names,
  c("DPL-1_Young-Adult", "EFL-1_Young-Adult", "PHA-4_Young-Adult",
    "W03F9.2_L4-Young-Adult-stage-larvae"))


# Gets the ChIP peaks upstream of all genes.
get.chip.locs = function(m) {

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

if (FALSE) {
output.dir = "git/sort_paper/tf/allResults/chip/"
system(paste("mkdir -p", output.dir))

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

if (TRUE) {
e = compute.enrichment.diff.cutoffs(chip.names[1:3], get.chip.locs,
  facs.enriched.depleted, c(1,2,3), c(0, 0.5, 0.7, 0.9), c(0))
e[,,,,,"p"][ is.na(e[,,,,,"p"]) ] = 1
e[,,,,,"p.corr"] = array(p.adjust(as.vector(e[,,,,,"p"]), method="fdr"),
    dim=dim(e[,,,,,"p.corr"]), dimnames=dimnames(e[,,,,,"p.corr"]))
enrich = e
save(enrich, file="git/sort_paper/tf/allResults/chip/facs.Rdata")

}


