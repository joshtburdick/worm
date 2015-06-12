# Basic stats about motifs.

load("git/sort_paper/tf/motif/hughes/motifCluster.Rdata")

# counts of Hughes motifs
# XXX this is hacky
count.hughes.motifs = function(org) {
  enrich = NULL
  load(paste0("git/sort_paper/tf/motif/hyperg/allResults2/",
    org, "_1.02/hier.300.clusters.Rdata"))
  m = unique(dimnames(enrich)[[1]])
  if (org != "Ce") {
    m = intersect(m, nr.motifs[[org]])
  }
  length(m)
}
hughes.motif.counts =
  sapply(c("Ce", "Dm", "Hs", "Mm"), count.hughes.motifs)
print(hughes.motif.counts)
cat("total Hughes motifs =", sum(hughes.motif.counts))

# ChIP signals
load("git/sort_paper/tf/motif/hyperg/allResults/chip/hier.300.clusters.Rdata")
chip.signals = dimnames(enrich)[[1]]
min.p.chip = apply(enrich[,,"p.corr",,,], 1, min)


if (FALSE) {
cat("number of 5kb motifs =", length(motifs.5kb), "\n")
cat("least significant 5kb =", max(min.p.5kb), "\n")

cat("number of Hughes motifs =", length(motifs.hughes), "\n")
cat("least significant Hughes =", max(min.p.hughes), "\n")

cat("total =", (length(motifs.5kb) + length(motifs.hughes)), "\n")
cat("total <= 1e-10 =", sum(c(min.p.5kb, min.p.hughes) <= 1e-10), "\n")

cat("number of ChIP signals =", length(chip.signals), "\n")
cat("least significant ChIP =", max(min.p.chip), "\n")
cat("num. ChIP significant =", sum(min.p.chip <= 0.05), "\n")
cat("num. ChIP < 1e-10 =", sum(min.p.chip <= 1e-10), "\n")
}





