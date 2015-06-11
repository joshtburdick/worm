# Shortens the list of enriched motifs.

library("MotIV")

# load("git/tf/motif/meme.format.pwm.Rdata")
load("git/sort_paper/tf/motif/hughes/motifCluster.Rdata")

# Use the existing motif clustering, and pick a motif per cluster
# (preferably a Ce motif.)
# XXX this may not be used
get.reduced.motifs = function() {
  a = cutree(hughes.motif.cluster[["all"]], h=0.2)
  a = a[ order( ! (names(a) %in% hughes.motif.cluster[["Ce"]]$labels), names(a)) ]

  # pick a motif per cluster
  m = a[ !duplicated(a) ]

  # convert from "group number" to "representative motif name"
  a1 = names(m)[ match(a, m) ]
  names(a1) = names(a)
  a1
}
motif.clustering = cutree(hughes.motif.cluster[["all"]], h=0.05)

# Clusters motifs, and then, for each cluster:
# - finds the most significant enrichment
# - combines all the genes which were called as orthologous, and
# - shows that list of genes (with NHRs collapsed, etc.)
# Args:
#   m - a table of motif enrichment
#   cl - the cluster number
motif.enrich.compact = function(m, cl) {

  # cluster the motifs, if there are multiple such
  if (nrow(m) > 1) {
    m$cluster = motif.clustering[ m$motif ]
  }
  # if there's only one motif, don't try to cluster it
  else {
    if (nrow(m) == 1)
      m$cluster = 1
  }
# browser()
  # find most significant enrichment for each cluster
  m = m[ order(m$p.corr) , ]
  m1 = m[ !duplicated(m$cluster) , ]

  if (nrow(m) > 0) {

    # add information about related orthologs
    m1$ortholog = ""
# browser()
    for(i in 1:nrow(m1)) {
      motifs.in.cl = m[ m[,"cluster"] == m1[i,"cluster"] , "motif" ]
      tf = unique(as.character(c(
        orthologs.by.motif[motifs.in.cl], recursive=TRUE)))
      m1[i,"ortholog"] = tf.list.annotate(tf, cl)
    }
  }
  m1
}

