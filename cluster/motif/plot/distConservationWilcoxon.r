# Compares motif location and conservation between genes
# in a cluster (and not), using a Wilcoxon test.

source("git/utils.r")

motif.gene.dir = "git/cluster/motif/distAndConservation/5kb/"

clusters = read.tsv(
  "git/cluster/hierarchical/hier.300.clusters/clusters.tsv")

# p-values
load("git/sort_paper/tf/motifEnrichment/hier.300.clusters.Rdata")

motif.enriched = read.tsv(gzfile(
  "git/sort_paper/tf/summary/motif/hier.300.clusters_merged.tsv.gz"))

# significance cutoff (since this includes all the enrichments)
motif.enriched = motif.enriched[ motif.enriched$p.corr <= 0.05 , ]

# sorting (to process most significant first)
motif.enriched = motif.enriched[ order(motif.enriched$p.corr) , ]

# data for histogram of amount of upstream conservation
upstream.cons.hist = read.tsv(gzfile(
  "git/tf/motif/conservation/cons_hist_WS220_5kb_upstream.tsv.gz"))

# Does a Mann-Whitney test for different distribution
# of upstream distance and conservation, between genes
# in the cluster, and outside the cluster.
compute.dist.conservation.wilcoxon = function(motif.enriched, output.file) {
  a = motif.enriched
  a$dist.U = NA
  a$dist.p = NA
  a$dist.p.corr = NA
  a$cons.U = NA
  a$cons.p = NA
  a$cons.p.corr = NA

  for(i in 1:2) {    # nrow(a)) {

    m = a[i,"motif"]
    cl = a[i,"group"]
    write.status(paste(i, m, cl))

    # read information for that motif
    r = read.table(paste(motif.gene.dir, m, "_upstreamMotifCons.tsv.gz", sep=""),
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
      a[i,"dist.p"] = w.dist$p.value

      w.cons = wilcox.test(r1$motif.cons, r0$motif.cons)
      a[i,"cons.U"] = w.cons$statistic
      a[i,"cons.p"] = w.cons$p.value
    }

    if (i %% 1000 == 0) {
      write.tsv(a, output.file)
    }
  }

  a$dist.p.corr = p.adjust(a$dist.p, method="fdr")
  a$cons.p.corr = p.adjust(a$cons.p, method="fdr")
  write.tsv(a, output.file)
}

if (TRUE) {
  compute.dist.conservation.wilcoxon(motif.enriched,
    gzfile("git/cluster/motif/plot/distConservationWilcoxon.tsv.gz"))
}

