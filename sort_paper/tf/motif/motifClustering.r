# Clusters motifs.

source("git/utils.r")

library("MotIV")

load("git/tf/motif/meme.format.pwm.Rdata")

# XXX hack to get motif names to use
r = read.tsv(gzfile("git/sort_paper/tf/motif/hyperg/table/hier.300.clusters.tsv.gz"))
m1 = unique(r$motif)

motifs = meme.format.pwm[m1]
dists = motifDistances(motifs)
motif.clusters = hclust(dists)

save(motif.clusters, file="git/sort_paper/tf/motif/motifClustering.Rdata")

