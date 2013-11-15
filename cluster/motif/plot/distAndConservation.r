# Plots how far away motifs are, and conservation of motifs,
# for various motifs enriched upstream of clusters.

source("git/utils.r")

output.dir = "git/cluster/motif/plot/distAndConservation/"

motif.gene.dir = "git/cluster/motif/distAndConservation/"

clusters = motifs.and.clusters = read.tsv(
  "git/cluster/hierarchical/hier.200.clusters/clusters.tsv")

motifs.and.clusters = read.table(
  "git/cluster/hierarchical/hier.200.clusters/uniqueKnownMotifEnrichment_5kb_0.5cons.tsv",
  header=TRUE, as.is=TRUE)

# Does a Mann-Whitney test for different distribution
# of upstream distance and conservation, between genes
# in the cluster, and outside the cluster.
compute.dist.conservation.wilcoxon = function() {
  r.dist = NULL
  r.conservation = NULL

  for(i in 1:5) {    # nrow(motifs.and.clusters)) {

    m = motifs.and.clusters[i,"motif"]
    cl = motifs.and.clusters[i,"cluster"]
    cat(backspace.string, i, m, cl)

    r = read.table(paste(motif.gene.dir, m, "_upstreamMotifCons.tsv.gz", sep=""),
      as.is=TRUE)
    colnames(r) = c("region.chr", "region.a", "region.b", "gene", "score",
      "strand", "motif", "motif.a", "motif.b", "motif.id", "motif.cons")
    r$upstream.dist = ifelse(r$strand=="+",
      r$motif.a - r$region.b, r$region.a - r$motif.a)
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

  r.dist$p.fdr = p.adjust(r.dist$p, method="fdr")
  r.conservation$p.fdr = p.adjust(r.conservation$p, method="fdr")

  list(r.dist = r.dist, r.conservation = r.conservation)
}

# Plots all the graphs.
plot.all = function(wilcox.results) {
  system(paste("mkdir -p", output.dir))

  for(i in 1:5) {  # nrow(motifs.and.clusters)) {

    m = motifs.and.clusters[i,"motif"]
    cl = motifs.and.clusters[i,"cluster"]

    r = read.table(paste(motif.gene.dir, m, "_upstreamMotifCons.tsv.gz", sep=""),
      as.is=TRUE)
    colnames(r) = c("region.chr", "region.a", "region.b", "gene", "score",
      "strand", "motif", "motif.a", "motif.b", "motif.id", "motif.cons")
    r$upstream.dist = ifelse(r$strand=="+",
      r$motif.a - r$region.b, r$region.a - r$motif.a)
    r$in.cluster = clusters[r$gene,"cluster"]==cl
    r = r[ !is.na(r$in.cluster) , ]

    if (sum(r$in.cluster > 0)) {
#      png(paste(output.dir, "/", m, " ", cl, ".png", sep=""),
#        width=1050, height=700)
      pdf(paste(output.dir, "/", m, " ", cl, ".pdf", sep=""),
        width=9, height=6)

      par(mfrow=c(2,3))
      r0 = r[ !r$in.cluster , ]
      r1 = r[ r$in.cluster , ]

      s = paste("Motif", m, "near genes in cluster", cl)
      plot(r1$upstream.dist, r1$motif.cons, col="#ff000080", pch=20, 
        xlim=c(-5000,0), ylim=c(0,1),
        main=s, xlab="Location relative to TSS", ylab="Conservation (PhastCons)")
      hist(r1$upstream.dist, col="#ff0000a0", breaks=30, xlim=c(-5000,0),
        main=s, xlab="Location relative to TSS")
      hist(r1$motif.cons, col="#ff0000a0", breaks=30, xlim=c(0, 1),
        main=s, xlab="Conservation (PhastCons)")

      s = paste("Motif", m, "near genes not in cluster", cl)
      plot(r0$upstream.dist, r0$motif.cons, col="#0000ff80", pch=20,
        xlim=c(-5000,0), ylim=c(0,1),
        main=s, xlab="Location relative to TSS", ylab="Conservation (PhastCons)")
      hist(r0$upstream.dist, col="#0000ffa0", breaks=30, xlim=c(-5000,0),
        main=s, xlab="Location relative to TSS")
      hist(r0$motif.cons, col="#0000ffa0", breaks=30, xlim=c(0, 1),
        main=s, xlab="Conservation (PhastCons)")

      dev.off()
    }
  }
}

if (TRUE) {
  r = compute.dist.conservation.wilcoxon()

  write.table(rbind(cbind(number="upstream.dist", r$r.dist),
    cbind(number="conservation", r$r.conservation)),
    file="git/cluster/motif/plot/distConservationWilcoxon.tsv",
    col.names=TRUE, row.names=FALSE)



# plot.all()
}
