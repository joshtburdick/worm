# Clusters motifs using MotIV.

library("MotIV")

source("git/tf/motif/motiv.utils.r")

# read in the de novo motifs
load("git/tf/motif/de.novo.motifs.Rdata")

# read in the known motifs
known.motif.meme.file = paste("data/tf/meme/motif_databases/",
    c("jolma2013.meme",
    "JASPAR_CORE_2009_insects.meme",
    "JASPAR_CORE_2009_nematodes.meme",
    "JASPAR_CORE_2009_vertebrates.meme"), sep="")

# read in known motifs
# XXX this is more complicated than presumably it needs to be
known.motifs = c(
  read.meme.file(known.motif.meme.file[[1]]),
  read.meme.file(known.motif.meme.file[[2]]),
  read.meme.file(known.motif.meme.file[[3]]),
  read.meme.file(known.motif.meme.file[[4]]))

#  read.meme.file("git/tf/motif/hughes_motif.meme"))

motifs = c(known.motifs, de.novo.motifs)

dists = motifDistances(motifs)
motif.clusters = hclust(dists)

# save(motif.clusters, known.motifs, file="git/tf/motif/clusterUsingMotIV_20141006.Rdata")

