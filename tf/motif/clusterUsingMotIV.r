# Clusters motifs using MotIV.

library("MotIV")

source("git/tf/motif/motiv.utils.r")

# read in the de novo motifs
# load("git/tf/motif/de.novo.motifs.Rdata")

# read in the known motifs
known.motif.meme.file = paste("data/tf/meme/motif_databases/",
    c("jolma2013.meme",
    "JASPAR_CORE_2009_insects.meme",
    "JASPAR_CORE_2009_nematodes.meme",
    "JASPAR_CORE_2009_vertebrates.meme"), sep="")

# add in Hughes motifs
load("git/tf/motif/meme.format.pwm.Rdata")
hughes.motif.names = grep("M...._1\\.01", names(meme.format.pwm))
hughes.motifs = meme.format.pwm[ hughes.motif.names ]

# read in known motifs
# XXX this is more complicated than presumably it needs to be
known.motifs = c(
  read.meme.file(known.motif.meme.file[[1]]),
  read.meme.file(known.motif.meme.file[[2]]),
  read.meme.file(known.motif.meme.file[[3]]),
  read.meme.file(known.motif.meme.file[[4]]),
  hughes.motifs)

motifs = known.motifs   # omitting de novo for now   c(known.motifs, de.novo.motifs)

# dists = motifDistances(motifs)
# motif.clusters = hclust(dists)

# save(motif.clusters, known.motifs, file="git/tf/motif/clusterUsingMotIV_20141006.Rdata")

