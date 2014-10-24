# Filters motifs to reduce redundancy, by clustering, and picking
# one motif from each cluster.

source("git/utils.r")

# the names of motifs which were counted
load("git/tf/motif/motifCount/motif.counts.Rdata")

# the motif clustering
load("git/tf/motif/clusterUsingMotIV.Rdata")

motif.cl = cutree(motif.clusters, h=0.05)

r = data.frame(motif.cl, priority=0, stringsAsFactors=FALSE)

# hack to preferentially choose a known motif (particularly
# one from the HT-selex experiment)
r$priority = rep(0, length(motif.cl))
# names(priority) = names(motif.cl)
r$priority[
  rownames(r) %in% names(known.motifs) ] = 1
r$priority[
  grep("(full|DBD)_\\d$", rownames(r)) ] = 2
r$priority[
  grep("(full|DBD)$", rownames(r)) ] = 3
r = r[ order(r$motif.cl, -r$priority) , ]

# for each cluster, pick a "canonical" name
r1 = r[ !duplicated(r$motif.cl) , ]

canonical.gene.name = rownames(r1)
names(canonical.gene.name) = r1$motif.cl
r$canonical.name =
  canonical.gene.name[ as.character(r$motif.cl) ]

r = r[ , c(2,3) ]

write.tsv(r, "git/tf/motif/motifFilter.tsv")

