# Plots how far away motifs are, and conservation of motifs,
# for various motifs enriched upstream of clusters.

source("git/utils.r")

output.dir = "git/cluster/motif/plot/distAndConservation/"

motif.gene.dir = "git/cluster/motif/distAndConservation/5kb/"

clusters = motifs.and.clusters = read.tsv(
  "git/cluster/hierarchical/hier.200.clusters/clusters.tsv")

motifs.and.clusters = read.table(
  "git/cluster/hierarchical/hier.200.clusters/uniqueKnownMotifEnrichment_5kb_0.5cons.tsv",
  header=TRUE, as.is=TRUE)

# data for histogram of amount of upstream conservation
upstream.cons.hist = read.tsv(gzfile(
  "git/tf/motif/conservation/upstreamRegionsWS220_conservationHist.txt.gz"))

# comparisons of upstream motif location and conservation
dist.cons.wilcoxon = read.table(
  "git/cluster/motif/plot/distConservationWilcoxon.tsv",
  sep="\t", header=TRUE, row.names=NULL, as.is=TRUE)

# Computes histogram counts "by hand". This avoids the issue
# that "hist(..., plot=FALSE)" is finicky about where the
# buckets are. Also, it's useful for making cumulative histograms.
hist.count = function(x, lim, num.buckets=20) {
  x = x[ !is.na(x) ]

  # scale to [0,1]
  x1 = floor( num.buckets * (x - lim[1]) / (lim[2] - lim[1]) )
  x1[ x1 < 0 ] = 0
  x1[ x1 >= num.buckets ] = num.buckets - 1

  counts = table(x1)

  r = rep(0, num.buckets)
  r[as.numeric(names(counts))+1] = as.vector(counts)

  data.frame(breaks = lim[1] + (c(0:(num.buckets-1)) / num.buckets) *
    (lim[2] - lim[1]),
    counts = r)
}

# Gets a histogram of amount of conservation.
# Args:
#   genes - genes to get conservation for
# Returns: data frame with columns
#   conservation - bucketized amount of conservation
#   bases - number of bases with that amount of conservation
get.conservation.buckets = function(genes) {
  g = intersect(genes, rownames(upstream.cons.hist))
  cons = apply(upstream.cons.hist[g,], 2, sum)

  data.frame(conservation = as.numeric(names(cons)),
    bases = as.vector(cons))
}

# Gets the number of base pairs upstream of each gene.
# (I'm using the conservation data here, which is admittedly
# sort of a hack.)
get.num.bp.upstream = function(genes) {
  g = intersect(genes, rownames(upstream.cons.hist))
  sizes = -apply(upstream.cons.hist[g,], 1, sum)

  h = hist.count(sizes, c(-5000,0), num.buckets=20)

  data.frame(bases = h$breaks, count = 250 * cumsum(h$counts))
}

# Plots all the graphs.
plot.all = function(wilcox.results) {
  system(paste("mkdir -p", output.dir))

  for(i in 1:nrow(motifs.and.clusters)) {

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
#        width=1050, height=1200)
      pdf(paste(output.dir, "/", m, " ", cl, ".pdf", sep=""),
        width=9, height=12)

      mat = matrix(c(1,1,1,2,3,4,5,6,7,8,9,10,11,12,13), nrow=5, byrow=TRUE)
      layout(mat, heights=c(0.5,1,1,1))

      # add labels
      par(mar=c(0,0,0,0) + 0.1)     # XXX
      plot(0,0, xlim=c(-1,1), ylim=c(0,5), type="n", bty="n",
        xaxt="n", yaxt="n", xlab="", ylab="")
      text(0,3, paste("Motif", m, "near genes in cluster", cl), cex=2.5)
#      text(0,0, "testing")
      text(0,1, paste("Motif enrichment p =",
        signif(motifs.and.clusters[i,"p.fdr"], 2)), cex=1.5)

      par(mar=c(5,4,4,2)+0.1)

      # par(mfrow=c(2,3))
      r0 = r[ !r$in.cluster , ]
      r1 = r[ r$in.cluster , ]

      # numbers for motifs near genes in clusters
      s = paste("Motif near genes in cluster")
      plot(r1$upstream.dist, r1$motif.cons, col="#ff000040", pch=20, 
        xlim=c(-5000,0), ylim=c(0,1),
        main=s, xlab="Location relative to TSS", ylab="Conservation (PhastCons)")
      hist(r1$upstream.dist, col="#ff0000a0", xlim=c(-5000,0),
        main=s, xlab="Location relative to TSS")
      dist.p = dist.cons.wilcoxon[ dist.cons.wilcoxon$motif==m &
        dist.cons.wilcoxon$cluster==cl &
        dist.cons.wilcoxon$number=="upstream.dist", "p.fdr" ]
      mtext(paste("Wilcoxon p =", signif(dist.p,2) ), cex=0.8)
      hist(r1$motif.cons, col="#ff0000a0", xlim=c(0, 1),
        main=s, xlab="Conservation (PhastCons)")
      cons.p = dist.cons.wilcoxon[ dist.cons.wilcoxon$motif==m &
        dist.cons.wilcoxon$cluster==cl &
        dist.cons.wilcoxon$number=="conservation", "p.fdr" ]
      mtext(paste("Wilcoxon p =", signif(cons.p,2)), cex=0.8)

      # similarly, numbers for motifs in other clusters
      s = paste("Motif near genes not in cluster")
      plot(r0$upstream.dist, r0$motif.cons, col="#0000ff40", pch=20,
        xlim=c(-5000,0), ylim=c(0,1),
        main=s, xlab="Location relative to TSS", ylab="Conservation (PhastCons)")
      hist(r0$upstream.dist, col="#0000ffa0", xlim=c(-5000,0),
        main=s, xlab="Location relative to TSS")
      hist(r0$motif.cons, col="#0000ffa0", xlim=c(0, 1),
        main=s, xlab="Conservation (PhastCons)")

      # same, for all sequence upstream of genes near clusters
      plot.new()

      bp = get.num.bp.upstream(clusters[ clusters[ r$gene,"cluster"]==cl , "gene" ])
      barplot(bp$count, col="#ff0000a0", space=0,
        main="All sequence near genes in cluster",
        xlab="Location relative to TSS", ylab="Frequency")
      # XXX
      par(new=TRUE)
      plot(0,0, type="n", xlim=c(-5000,0), ylim=c(0,1), bty="n",
        xaxt="n", yaxt="n", xlab="", ylab="")
      axis(1, line=0.4)

      cons = get.conservation.buckets(
        clusters[ clusters[ r$gene,"cluster"]==cl , "gene" ])
      barplot(cons$bases, col="#ff0000a0", space=0,
        main="All sequence near genes in cluster",
        xlab="Conservation (PhastCons)", ylab="Frequency")
      # XXX
      par(new=TRUE)
      plot(0,0, type="n", xlim=c(0,1), ylim=c(0,1), bty="n",
        xaxt="n", yaxt="n", xlab="", ylab="")
      axis(1, line=0.4)


      # and again, background numbers in all clusters
      plot.new()

      bp = get.num.bp.upstream(clusters[ clusters[ r$gene,"cluster"]!=cl , "gene" ])
      barplot(bp$count, col="#0000ffa0", space=0,
        main="All sequence near genes not in cluster",
        xlab="Location relative to TSS", ylab="Frequency")
      # XXX
      par(new=TRUE)
      plot(0,0, type="n", xlim=c(-5000,0), ylim=c(0,1), bty="n",
        xaxt="n", yaxt="n", xlab="", ylab="")
      axis(1, line=0.4)

      cons = get.conservation.buckets(
        clusters[ clusters[ r$gene,"cluster"]!=cl , "gene" ])
      barplot(cons$bases, col="#0000ffa0", space=0,
        main="All sequence near genes not in cluster",
        xlab="Conservation (PhastCons)", ylab="Frequency")
      # XXX
      par(new=TRUE)
      plot(0,0, type="n", xlim=c(0,1), ylim=c(0,1), bty="n",
        xaxt="n", yaxt="n", xlab="", ylab="")
      axis(1, line=0.4)

      dev.off()
    }
  }
}

if (TRUE) {

  plot.all()
}

