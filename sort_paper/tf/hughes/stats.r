# Stats about Hughes motifs.

source("git/utils.r")
load("git/sort_paper/tf/motif/hughes/motifMatrix.Rdata")

motifs.with.pwms =
  unique(c(sapply(hughes.motif.matrix, names), recursive=TRUE))

motifs.with.counts = sub(".Rdata$", "",
  list.files("git/sort_paper/tf/motif/upstreamCount/Ce_1.02/"))


print.stats = function(f) {
  a = read.table(paste0("data/tf/hughes/Ce_1.02/", f),
    sep="\t", header=TRUE, as.is=TRUE)
  a = a[ a$Motif_ID != "." , ]
  cat("------", f, "\n")
  cat("num genes with motifs = ", length(unique(a$TF_Name)), "\n")
  cat("num motifs = ", length(unique(a$Motif_ID)), "\n")

  cat("[after filtering out motifs w/o counts (TRANSFAC, or v. degenerate)]\n")
  a = a[ a$Motif_ID %in% motifs.with.counts , ]  
  cat("num genes with motifs = ", length(unique(a$TF_Name)), "\n")
  cat("num motifs = ", length(unique(a$Motif_ID)), "\n\n")
}

print.stats("TF_Information.txt")
print.stats("TF_Information_all_motifs.txt")
print.stats("TF_Information_all_motifs_plus.txt")

