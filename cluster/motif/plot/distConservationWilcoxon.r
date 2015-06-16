# Compares motif location and conservation between genes
# in a cluster (and not), using a Wilcoxon test.

source("git/utils.r")

motif.gene.dir = "/media/jburdick/disk2/jburdick/distAndConservation/"

clusters = read.tsv(
  "git/cluster/hierarchical/hier.300.clusters/clusters.tsv")

# p-values
# motif.enriched = read.tsv(gzfile(
#   "git/sort_paper/tf/summary/motif/hier.300.clusters_merged.tsv.gz"))
motif.enriched = read.tsv(gzfile(
  "git/sort_paper/tf/motif/hyperg/summary/hughes/hier.300.clusters.tsv.gz"))

# sorting (to process most significant first)
motif.enriched = motif.enriched[ order(motif.enriched$p.corr) , ]

# only including most significant enrichments
motif.enriched = motif.enriched[1:10000,]

# data for histogram of amount of upstream conservation
# upstream.cons.hist = read.tsv(gzfile(
#   "git/tf/motif/conservation/cons_hist_WS220_5kb_upstream.tsv.gz"))

# XXX hack to get the appropriate motif from whichever directory
get.motif.filename = function(m) {
  for(org in c("Ce", "Dm", "Hs", "Mm")) {
    f = paste0(motif.gene.dir, org, "_1.02/", m,
      "_upstreamMotifCons.tsv.gz")
    if (file.exists(f))
      return(f)
  }

  cat("\nfailed to find info for motif", m, "\n")
  return(NA)
}

# Does a Mann-Whitney test for different distribution
# of upstream distance and conservation, between genes
# in the cluster, and outside the cluster.
compute.dist.conservation.wilcoxon = function(motif.enriched, output.file) {
  a = motif.enriched
  a$dist.U = NA
  a$dist.closer = NA
  a$dist.p = NA
  a$dist.p.corr = NA
  a$cons.U = NA
  a$cons.higher = NA
  a$cons.p = NA
  a$cons.p.corr = NA

  for(i in 1:nrow(a)) {

    m = a[i,"motif"]
    cl = a[i,"group"]
    write.status(paste(i, m, cl))

    # read information for that motif
    r = read.table(get.motif.filename(m),
      as.is=TRUE)
    colnames(r) = c("region.chr", "region.a", "region.b", "gene", "score",
      "strand", "motif", "motif.a", "motif.b", "motif.id", "motif.cons")
    r$upstream.dist = ifelse(r$strand=="+",
      r$motif.a - r$region.b, r$region.a - r$motif.a)
    r$upstream.dist[ r$upstream.dist > 0 ] = 0

    # compute genes in this cluster
    r$in.cluster = clusters[r$gene,"cluster"]==cl
    r = r[ !is.na(r$in.cluster) , ]

    if (sum(r$in.cluster) > 0 && sum(!r$in.cluster) > 0) {
      r0 = r[ !r$in.cluster , ]
      r1 = r[ r$in.cluster , ]

      w.dist = wilcox.test(r1$upstream.dist, r0$upstream.dist)
      a[i,"dist.U"] = w.dist$statistic
      a[i,"dist.closer"] = w.dist$statistic > 0.5 * nrow(r0) * nrow(r1)
      a[i,"dist.p"] = w.dist$p.value

      w.cons = wilcox.test(r1$motif.cons, r0$motif.cons)
      a[i,"cons.U"] = w.cons$statistic
      a[i,"cons.higher"] = w.cons$statistic > 0.5 * nrow(r0) * nrow(r1)
      a[i,"cons.p"] = w.cons$p.value
    }

    if (i %% 100 == 0) {
      write.tsv(a, gzfile(output.file))
    }
  }

  a$dist.p.corr = p.adjust(a$dist.p, method="fdr")
  a$cons.p.corr = p.adjust(a$cons.p, method="fdr")
  write.tsv(a, gzfile(output.file))
}

if (TRUE) {
  compute.dist.conservation.wilcoxon(motif.enriched,
    "git/cluster/motif/plot/distConservationWilcoxon.tsv.gz")
}

