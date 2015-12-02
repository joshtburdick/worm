# Plots how far away motifs are, and conservation of motifs,
# for various motifs enriched upstream of clusters.

source("git/utils.r")
source("git/plot/utils.r")
source("git/sort_paper/tf/motif/hughes/motifInfo.r")

output.dir = "git/cluster/motif/plot/distAndConservation2/"

# motif.gene.dir = "git/cluster/motif/distAndConservation/5kb/"
motif.gene.dir = "/media/jburdick/data1/jburdick/distAndConservation/"

# this will define which genes are in the cluster
clusters = read.tsv(
  "git/cluster/hierarchical/hier.300.clusters/clusters.tsv")

# enrichment of each motif
motif.enrich = read.tsv(gzfile(
  "git/sort_paper/tf/summary/motif/hier.300.clusters_merged.tsv.gz"))

# data for histogram of amount of upstream conservation
upstream.cons.hist = read.tsv(gzfile(
  "git/tf/motif/conservation/cons_hist_WS220_5kb_upstream.tsv.gz"))

# comparisons of upstream motif location and conservation
dist.cons.wilcoxon = read.tsv(gzfile(
  "git/cluster/motif/plot/distConservationWilcoxon.tsv.gz"))
dist.cons.wilcoxon$dist.p.corr =
  p.adjust(dist.cons.wilcoxon$dist.p, method="fdr")
dist.cons.wilcoxon$cons.p.corr =
  p.adjust(dist.cons.wilcoxon$cons.p, method="fdr")

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

# get amount of upstream sequence for each gene
upstream.dist = {
  upstream.bed = read.table("git/tf/motif/motifCount/regions/upstream_5kb_WS220.bed", as.is=TRUE)
  colnames(upstream.bed) = c("chr", "a", "b", "gene", "score", "strand")
  upstream.bed$upstream.dist = upstream.bed$b - upstream.bed$a
  upstream.dist = upstream.bed$upstream.dist
  names(upstream.dist) = upstream.bed$gene
  upstream.dist
}

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

# Plots distribution of motifs upstream of genes.
# Args: a, one-row data.frame including fields:
#   motif - the motif
#   group - the cluster number
#   p.corr - p-value for enrichment of this motif-cluster pair
#   dist.p.corr - p-value for motifs in cluster being nearer TSS
#   cons.p.corr - p-value for motifs in cluster being more conserved
# (all p-values are presumed to have been corrected already)
# Side effects: plots the graph.
plot.motif.loc.dist = function(a, motif.score.cutoff=40) {

  m = a$motif
  cl = as.character(a$group)
  enrich.p = a$p.corr

  system(paste("mkdir -p", output.dir))

  r = read.table(get.motif.filename(m), as.is = TRUE)

  colnames(r) = c("region.chr", "region.a", "region.b", "gene", "score",
    "strand", "motif.chr", "motif.a", "motif.b", "motif.id",
"motif.score", "motif.strand", "motif.cons")
  r$upstream.dist = ifelse(r$strand=="+",
    r$motif.a - r$region.b, r$region.a - r$motif.a)
  # XXX this is presumably due to really long motifs near
  # where a gene starts
  r$upstream.dist[ r$upstream.dist > 0 ] = 0

  # tack on normalized upstream location
  r$upstream.dist.normalized =
    r$upstream.dist / upstream.dist[ r$gene ]
  r$upstream.dist.normalized[ r$upstream.dist.normalized < -1 ] = -1
  r$upstream.dist.normalized[ r$upstream.dist.normalized > 0 ] = 0
# browser()

  r$in.cluster = clusters[r$gene,"cluster"]==cl
  r = r[ !is.na(r$in.cluster) , ]

  # XXX this can also affect these graphs
  r = r[ r$motif.score >= motif.score.cutoff , ]

  if (sum(r$in.cluster > 0)) {
#      png(paste(output.dir, "/", m, " ", cl, ".png", sep=""),
#        width=1050, height=1200)
    pdf(paste(output.dir, "/", m, " ", cl, ".pdf", sep=""),
      width=9, height=6)    # was "height=12"

#    mat = matrix(c(1,1,1,2,3,4,5,6,7,8,9,10,11,12,13), nrow=5, byrow=TRUE)
#    layout(mat, heights=c(0.5,1,1,1))
    mat = matrix(c(1,1,1,1, 2,3,4,5, 6,7,8,9), nrow=3, byrow=TRUE)
    layout(mat, heights=c(0.2,1,1), widths=c(0.15,1,1,1))

    # add labels
    par(mar=c(0,0,0,0) + 0.1)     # XXX
    plot(0,0, xlim=c(-1,1), ylim=c(0,5), type="n", bty="n",
      xaxt="n", yaxt="n", xlab="", ylab="")

print(motif.info[m,])
    motif.name.1 = as.expression(motif.info[m, "motif.name"])

    # XXX hacks
    if (motif.info[m,"motif.name"] == "RFX3_1")
      motif.name.1 = expression(italic("daf-19"))
    if (m == "M2336_1.02")
      motif.name.1 = expression(italic("eor-1"))
    if (m == "M0892_1.02")
      motif.name.1 = expression(italic("ceh-2"))

#    if (m == "M1906_1.02")
#      motif.name.1 = expression(italic("klf-2"))

    text(0,3,
      expr.format(expression(m * " motif near genes in cluster " * cl),
        list(m=motif.name.1, cl=cl)), cex=2)
 
#    text(0,3, paste(motif.name.1, "motif", "near genes in cluster", cl), cex=2)

#    text(0,1, paste("Motif enrichment p =",
#      signif(enrich.p, 2)), cex=1.5)
    text(0,1, expr.format(expression("Motif enrichment p " * p),
      list(p = format.p(signif(enrich.p,2), include.equals.sign=TRUE))), cex=1.5)

    par(mar=c(5,4,4,1)+0.1)

    # par(mfrow=c(2,3))
    r0 = r[ !r$in.cluster , ]
    r1 = r[ r$in.cluster , ]

    par(mar=c(0,0,0,0))
    plot.new()

    s = "" # paste("Motif near genes in cluster")

    # numbers for motifs near genes in clusters
if (FALSE) {

    plot(r1$upstream.dist, r1$motif.cons, col="#ff000040",
      xlim=c(-5000,0), ylim=c(0,1),
      pch=183, font=5, xaxt="n", yaxt="n",
      main=s, xlab="Location relative to TSS", ylab="Conservation (PhastCons)")
    axis(1)
    axis(2)
}
    par(mar=c(5,4,2,1)+0.1, cex.lab=1.2)
    hist(r1$upstream.dist.normalized, col="#ff0000a0", xlim=c(-1,0),
      main=s, xlab="")
    mtext(expression(bold("Genes in cluster")),
      side=2, line=5, cex=1, xpd=NA)

    hist(r1$upstream.dist, col="#ff0000a0", xlim=c(-5000,0),
      main=s, xlab="")

#    dist.p = dist.cons.wilcoxon[ dist.cons.wilcoxon$motif==m &
#      dist.cons.wilcoxon$cluster==cl &
#      dist.cons.wilcoxon$number=="upstream.dist", "p.fdr" ]
    mtext(expr.format(expression("Wilcoxon p " * p),
      list(p = format.p(signif(a$dist.p.corr,2), include.equals.sign=TRUE))), cex=0.8)

    hist(r1$motif.cons, col="#ff0000a0", xlim=c(0, 1),
      main=s, xlab="")
#    cons.p = dist.cons.wilcoxon[ dist.cons.wilcoxon$motif==m &
#      dist.cons.wilcoxon$cluster==cl &
#      dist.cons.wilcoxon$number=="conservation", "p.fdr" ]
    mtext(expr.format(expression("Wilcoxon p " * p),
      list(p = format.p(signif(a$cons.p.corr,2), include.equals.sign=TRUE))), cex=0.8)

    par(mar=c(0,0,0,0))
    plot.new()

    # similarly, numbers for motifs in other clusters
    s = "" # paste("Motif near genes not in cluster")
if (FALSE) {
    plot(r0$upstream.dist, r0$motif.cons, col="#0000ff40",
      pch=183, font=5, xaxt="n", yaxt="n",
      xlim=c(-5000,0), ylim=c(0,1),
      main=s, xlab="Location relative to TSS", ylab="Conservation (PhastCons)")
    axis(1)
    axis(2)
}
    par(mar=c(5,4,2,1)+0.1, cex.lab=1.2)
    hist(r0$upstream.dist.normalized, col="#0000ffa0", xlim=c(-1,0),
      main=s, xlab=expression(bold("Normalized location relative to 5' end")))
    mtext(expression(bold("Genes not in cluster")),
      side=2, line=5, cex=1, xpd=NA)

    hist(r0$upstream.dist, col="#0000ffa0", xlim=c(-5000,0),
      main=s, xlab=expression(bold("Location relative to transcript 5' end")))

    hist(r0$motif.cons, col="#0000ffa0", xlim=c(0, 1),
      main=s, xlab=expression(bold("Conservation (PhastCons)")))

    dev.off()
  }

}

# Plots all the graphs (deprecated.)
plot.all.old = function(wilcox.results) {
  system(paste("mkdir -p", output.dir))

  for(i in 1:nrow(motifs.and.clusters)) {

    m = motifs.and.clusters[i,"motif"]
    cl = motifs.and.clusters[i,"cluster"]

# FIXME
# plot.motif.loc.dist(...)

  }
}

# Plots one of these.
# Args:
#   motif - the motif
#   cluster - the cluster
# Side effects: plots the graphs, in the output directory.
plot.one = function(motif, cluster, motif.score.cutoff=40) {
  a = dist.cons.wilcoxon[dist.cons.wilcoxon$motif == motif &
    dist.cons.wilcoxon$group == cluster , ]
  stopifnot(nrow(a)==1)
  plot.motif.loc.dist(a, motif.score.cutoff)
}

plot.some = function() {
  plot.one("RFX2_DBD", "286")
  plot.one("M0665_1.00", "30")
  plot.one("M1936_1.00", "30", motif.score.cutoff=30)

#  for(i in 1:50) {
#    plot.motif.loc.dist(dist.cons.wilcoxon[i,], motif.score.cutoff=40)
#  }
}

# Plots cases of motifs which seem to be distant (but conserved).
plot.conserved.distant = function() {
  a = dist.cons.wilcoxon
  distant.motifs = a[ ( (a$p.corr <= 1e-10) &
    (!a$dist.closer) & (a$dist.p.corr <= 1e-5) &
    (a$cons.higher) & (a$cons.p.corr <= 1e-2) ) , ]
# distant.motifs = distant.motifs[ !duplicated(distant.motifs$group) , ]

#  browser()
  print(dim(distant.motifs))
  for(i in 1:nrow(distant.motifs)) {
    cat(i, "\n")
    print(distant.motifs[i,])
    plot.one(distant.motifs[i,"motif"], distant.motifs[i,"group"])
  }
  distant.motifs
}

if (TRUE) {
  plot.one("M2331_1.02", "284")
  plot.one("M0892_1.02", "284")

  # conserved, and further away
  plot.one("M2336_1.02", "284")
#  plot.one("M4567_1.02", "46")
#  plot.one("M1906_1.02", "245")

  # conserved, and closer
  plot.one("M5775_1.02", "286")
#  plot.some()
}

# distant.motifs = plot.conserved.distant()

