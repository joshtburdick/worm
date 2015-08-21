# Looks for pairs of motifs which are both enriched
# upstream of some cluster.

source("git/utils.r")
source("git/tf/motif/enrichment/motifHyperg.r")

# a clustering
cl = {
  cl1 = read.tsv(
    "git/cluster/hierarchical/hier.300.clusters/clusters.tsv")
  cl = cl1$cluster
  names(cl) = cl1$gene
  cl
}

# motif.gene.dir = "git/cluster/motif/distAndConservation/5kb/"
motif.gene.dir = "/media/jburdick/disk2/jburdick/distAndConservation/"





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

# table of which cutoffs were significant the most often
if (TRUE) {
mh = read.tsv(
  "git/sort_paper/tf/motif/hyperg/summary/hughes/hier.300.clusters.tsv.gz")
mh = mh[ mh$p.corr <= 0.05 , ]
m = unique(mh$motif)
most.frequent =
  function(x) as.character(names(sort(table(x), decreasing=TRUE)))[1]
most.significant.cutoffs = data.frame(
  upstream.dist.kb =
    c(by(mh$upstream.dist.kb, mh$motif, most.frequent))[m],
  conservation =
    c(by(mh$conservation, mh$motif, most.frequent))[m],
  motif.score =
    c(by(mh$motif.score, mh$motif, most.frequent))[m],
  stringsAsFactors=FALSE
)
}

# Gets upstream motif occurrences, at the thresholds which
# were most often most significant.
# Args:
#   m - name of the motif
# Returns: table of the upstream occurrences of that motif,
#   filtered by the above cutoffs.
get.motif.upstream.occurrences = function(m) {
  r = read.table(get.motif.filename(m), as.is = TRUE)
  colnames(r) =
    c("region.chr", "region.a", "region.b",
    "gene", "score", "strand",
    "motif", "motif.a", "motif.b", "motif.id", "motif.score",
    "motif.strand", "motif.cons")
  r$upstream.dist = ifelse(r$strand=="+",
    r$motif.b - r$region.b, r$region.a - r$motif.a)

  cutoffs = as.numeric(most.significant.cutoffs[m,])
  r = r[ (r$gene %in% names(cl) &
    r$upstream.dist >= -1000 * cutoffs[1]) &
    (r$motif.cons >= cutoffs[2]) &
    (r$motif.score >= cutoffs[3]) , ]

  r
}

# Finds the distance between the nearest occurrrences of
# two motifs.
# Args:
#   dist.1, dist.2 - the upstream distances for each motif
# Returns: vector of distances (for genes with both motifs)
motif.nearest.dist = function(dist.1, dist.2) {
  g = intersect(names(dist.1), names(dist.2))
  f = function(x)
    min(abs(outer(dist.1[[x]], dist.2[[x]], "-")))
  sapply(g, f)
}

# Checks for enrichment of a pair of motifs, and whether
# they're closer together than you'd expect by chance.
# Args:
#   m1, m2 - the motifs to check
# Returns: data frame of results for all clusters
pairwise.motif.enrich = function(m1, m2, dist.fun) {
  # get motif occurrences
  r1 = get.motif.upstream.occurrences(m1)
  r2 = get.motif.upstream.occurrences(m2)

  # count cluster sizes, and number of motif occurrences for each
  cluster.sizes = tapply(cl, cl, length)
  has.both.motifs = names(cl) %in% intersect(r1$gene, r2$gene)
  genes.with.motifs.in.cluster = tapply(has.both.motifs, cl, sum)

  # enrichment statistics
  enrich = (genes.with.motifs.in.cluster / cluster.sizes) /
    (sum(genes.with.motifs.in.cluster) / sum(cluster.sizes))
  p = motif.enrich.hyperg(genes.with.motifs.in.cluster, cluster.sizes,
    sum(genes.with.motifs.in.cluster), sum(cluster.sizes))

  # find motif locations
  ml1 = tapply(r1$upstream.dist, r1$gene, c)
  ml2 = tapply(r2$upstream.dist, r2$gene, c)

  motif.dist = motif.nearest.dist(ml1, ml2)
  cl1 = cl[ names(motif.dist) ]
  # ??? pad this with an arbitrary large number, for genes
  # without the pair of motifs?

  r = NULL
  clusters = sort(unique(cl))
  for(i in clusters) {
    write.status(paste(m1, m2, i))
    dist0 = motif.dist[ cl1 != i ]
    dist1 = motif.dist[ cl1 == i ]
    dist0.mean = NA
    dist1.mean = NA
    motif.dist.p = 1
    if ((length(dist0) > 0) && (length(dist1) > 0)) {
      dist0.mean = mean(dist0)
      dist1.mean = mean(dist1)
      w = wilcox.test(dist1, dist0, alternative="less")
      motif.dist.p = w$p.value
    }
    r = rbind(r, data.frame(
      m1 = m1, m2 = m2, cl = i,
      enrich = enrich[i], enrich.p = p[i],
      dist1 = dist1.mean, dist0 = dist0.mean,
      motif.dist.p = motif.dist.p, stringsAsFactors=FALSE))
  }

#  rownames(r) = cl
  r
}

# example, using two random motifs
# z = pairwise.motif.enrich("M0103_1.02", "M1865_1.02")

