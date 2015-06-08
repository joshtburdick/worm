# Shortens the list of enriched motifs.

library("MotIV")

load("git/tf/motif/meme.format.pwm.Rdata")

# Clusters motifs, and then, for each cluster:
# - finds the most significant enrichment
# - combines all the genes which were called as orthologous, and
# - shows that list of genes (with NHRs collapsed, etc.)
motif.enrich.compact = function(m, cl, dist.cutoff = 0.05) {

  # cluster the motifs, if there are multiple such
  if (nrow(m) > 1) {
    motifs = meme.format.pwm[ m$motif ]

    d = motifDistances(motifs)
    h = hclust(d, method="average")
    m$cluster = cutree(h, h = dist.cutoff)
  }
  # if there's only one motif, don't try to cluster it
  else {
    if (nrow(m) ==1)
      m$cluster = 1
  }

  # find most significant enrichment for each cluster
  m = m[ order(m$p.corr) , ]
  m1 = m[ !duplicated(m$cluster) , ]

  if (nrow(m) > 0) {

    # add information about related orthologs
    m1$ortholog = ""

    for(i in 1:nrow(m1)) {
  # browser()
      motifs.in.cl = m[ m[,"cluster"] == m1[i,"cluster"] , "motif" ]
      tf = unique(as.character(c(
        orthologs.by.motif[motifs.in.cl], recursive=TRUE)))
      m1[i,"ortholog"] = tf.list.annotate(tf, cl)
    }
  }
  m1
}

