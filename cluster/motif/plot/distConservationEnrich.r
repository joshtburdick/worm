# Compares motif location and conservation between genes
# in a cluster (and not), based on enrichment/depletion
# of motifs near the 5' transcript end, or in conserved
# regions.

source("git/utils.r")

motif.gene.dir = "/media/jburdick/disk2/jburdick/distAndConservation/"

clusters = read.tsv(
  "git/cluster/hierarchical/hier.300.clusters/clusters.tsv")

# p-values
# motif.enriched = read.tsv(gzfile(
#   "git/sort_paper/tf/summary/motif/hier.300.clusters_merged.tsv.gz"))
motif.enriched = read.tsv(gzfile(
  "git/sort_paper/tf/motif/hyperg/summary/hughes/hier.300.clusters.tsv.gz"))
motif.enriched = motif.enriched[ motif.enriched$p.corr <= 0.05 , ]

# sorting (to process most significant first)
motif.enriched = motif.enriched[ order(motif.enriched$p.corr) , ]

# only including most significant enrichments
# motif.enriched = motif.enriched[1:10000,]

# filter redundant motifs
load("git/sort_paper/tf/motif/hughes/motifCluster.Rdata")
motif.clustering = cutree(hughes.motif.cluster[["all"]], h=0.01)
motif.enriched = motif.enriched[
  ! duplicated(paste(motif.clustering[ motif.enriched$motif ],
   motif.enriched$group)) , ]

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

# Looks for enrichment of motif distance or conservation.
compute.dist.conservation.proportion = function(motif.enriched, output.file) {
  # cutoffs
  dist.cutoff = -1000
  cons.cutoff = 0.5

  a = motif.enriched

  a$num.motifs = NA

  a$dist.cluster = NA
  a$dist.enrich = NA
  a$dist.chisq = NA
  a$dist.p = NA
  a$dist.p.corr = NA

  a$cons.cluster = NA
  a$cons.enrich = NA
  a$cons.chisq = NA
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
      "strand", "motif",
      "motif.a", "motif.b", "motif.id", "motif.score", "motif.strand", "motif.cons")
    r$upstream.dist = ifelse(r$strand=="+",
      r$motif.a - r$region.b, r$region.a - r$motif.a)
    r$upstream.dist[ r$upstream.dist > 0 ] = 0
    r$in.cluster = clusters[r$gene,"cluster"]==cl
    r = r[ !is.na(r$in.cluster) , ]

    # count motifs in each category
    num.motifs = nrow(r)
    motifs.in.cluster = sum(r$in.cluster)

    if ((motifs.in.cluster > 0) && (motifs.in.cluster < num.motifs)) {
      a[i,"num.motifs"] = num.motifs

      # XXX these two are pretty repetitive
      closer = r[,"upstream.dist"] >= dist.cutoff
      dist.cluster = sum( closer[r$in.cluster] )
      dist.total = sum( closer )
      a[i,"dist.cluster"] = dist.cluster
      a[i,"dist.enrich"] = (dist.cluster / motifs.in.cluster) /
        (dist.total / num.motifs)
      t1 = chisq.test(c(dist.cluster, motifs.in.cluster),
        p=c(dist.total, num.motifs), rescale.p = TRUE)
      a[i,"dist.chisq"] = as.numeric(t1$statistic)
      a[i,"dist.p"] = t1$p.value

      cons = r[,"motif.cons"] >= cons.cutoff
      cons.cluster = sum( cons[r$in.cluster] )
      cons.total = sum( cons )
      a[i,"cons.cluster"] = cons.cluster
      a[i,"cons.enrich"] = (cons.cluster / motifs.in.cluster) /
        (cons.total / num.motifs)
      t2 = chisq.test(c(cons.cluster, motifs.in.cluster),
        p=c(cons.total, num.motifs), rescale.p = TRUE)
      a[i,"cons.chisq"] = as.numeric(t2$statistic)
      a[i,"cons.p"] = t2$p.value
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
  compute.dist.conservation.proportion(motif.enriched,
    "git/cluster/motif/plot/distConservationEnrich.tsv.gz")
}

