# Basic stats about motifs.

# original motifs
load("git/sort_paper/tf/motif/hyperg/allResults/5kb/hier.300.clusters.Rdata")
motifs.5kb = dimnames(enrich)[[1]]
min.p.5kb = apply(enrich[,,"p.corr",,,], 1, min)

# motifs from Hughes et al
load("git/sort_paper/tf/motif/hyperg/allResults/hughes/hier.300.clusters.Rdata")
motifs.hughes = dimnames(enrich)[[1]]
min.p.hughes = apply(enrich[,,"p.corr",,,], 1, min)

# ChIP signals
load("git/sort_paper/tf/motif/hyperg/allResults/chip/hier.300.clusters.Rdata")
chip.signals = dimnames(enrich)[[1]]
min.p.chip = apply(enrich[,,"p.corr",,,], 1, min)



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




