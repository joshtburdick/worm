# Compares motif location and conservation between genes
# in a cluster (and not), using a Wilcoxon test.

source("git/utils.r")

motif.gene.dir = "git/cluster/motif/distAndConservation/5kb/"

clusters = read.tsv(
  "git/cluster/hierarchical/hier.300.clusters/clusters.tsv")

# p-values
load("git/sort_paper/tf/motifEnrichment/hier.300.clusters.Rdata")

# data for histogram of amount of upstream conservation
upstream.cons.hist = read.tsv(gzfile(
  "git/tf/motif/conservation/cons_hist_WS220_5kb_upstream.tsv.gz"))

# Does a Mann-Whitney test for different distribution
# of upstream distance and conservation, between genes
# in the cluster, and outside the cluster.
compute.dist.conservation.wilcoxon = function() {
  r.dist = NULL
  r.conservation = NULL

  # find most significant p-value
  min.p = apply(enrich[,,,,,"p.corr"], c(1,2), min)

  for(i in 1:nrow(min.p)) {
    m = rownames(min.p)[i]

    # read information for that motif
    r = read.table(paste(motif.gene.dir, m, "_upstreamMotifCons.tsv.gz", sep=""),
      as.is=TRUE)
    colnames(r) = c("region.chr", "region.a", "region.b", "gene", "score",
      "strand", "motif", "motif.a", "motif.b", "motif.id", "motif.cons")
    r$upstream.dist = ifelse(r$strand=="+",
      r$motif.a - r$region.b, r$region.a - r$motif.a)
    r$upstream.dist[ r$upstream.dist > 0 ] = 0

    # loop through clusters that motif was associated with
    for(cl in colnames(min.p)[ which(min.p[i,] <= 0.05) ]) {
#      cl = motifs.and.clusters[i,"cluster"]
      write.status(paste(i, m, cl))

      # compute genes in this cluster
      r$in.cluster = clusters[r$gene,"cluster"]==cl
      r = r[ !is.na(r$in.cluster) , ]

      if (sum(r$in.cluster) > 0 && sum(!r$in.cluster) > 0) {
        r0 = r[ !r$in.cluster , ]
        r1 = r[ r$in.cluster , ]

        a = wilcox.test(r1$upstream.dist, r0$upstream.dist)
        a1 = data.frame(motif = m, cluster = cl,
          num.motifs.in.cluster = nrow(r1),
          num.motifs.in.background = nrow(r0),
          U = a$statistic, p = a$p.value, stringsAsFactors=FALSE)
        r.dist = rbind(r.dist, a1)

        a = wilcox.test(r1$motif.cons, r0$motif.cons)
        a1 = data.frame(motif = m, cluster = cl,
          num.motifs.in.cluster = nrow(r1),
          num.motifs.in.background = nrow(r0),
          U = a$statistic, p = a$p.value, stringsAsFactors=FALSE)
        r.conservation = rbind(r.conservation, a1)
      }
    }
  }

  r.dist$p.fdr = p.adjust(r.dist$p, method="fdr")
  r.conservation$p.fdr = p.adjust(r.conservation$p, method="fdr")

  list(r.dist = r.dist, r.conservation = r.conservation)
}

if (TRUE) {
r = compute.dist.conservation.wilcoxon()

write.table(rbind(cbind(number="upstream.dist", r$r.dist),
  cbind(number="conservation", r$r.conservation)),
  file="git/cluster/motif/plot/distConservationWilcoxon.tsv",
  sep="\t", col.names=TRUE, row.names=FALSE)
}

