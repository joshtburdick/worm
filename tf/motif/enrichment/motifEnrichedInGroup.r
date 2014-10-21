# Looks for motifs enriched in a group of genes.

source("git/utils.r")

load("git/tf/motif/count/upstreamRegionSize.Rdata")

stat.names = c("motifs.cluster", "motifs.background",
  "bp.cluster", "bp.background", "enrich", "chisq", "p", "p.corr")

# Computes enrichment for one motif, based on the relevant counts.
compute.enrichment = function(motifs.cluster, motifs.background,
    bp.cluster, bp.background) {

  # compute counts
  r = c(motifs.cluster = motifs.cluster,
    motifs.background = motifs.background,
    bp.cluster = bp.cluster,
    bp.background = bp.background,
    enrichment =
      (motifs.cluster / motifs.background) / (bp.cluster / bp.background),
    chisq = NA, p = 1, p.corr = NA)

  # only include enrichments
  if (!is.na(r["enrichment"]) && r["enrichment"] > 1) {

    # do the chi-square test; on error, just set the p-value to 1
    tryCatch({
      a = chisq.test(c(motifs.cluster, motifs.background),
        p = c(bp.cluster, bp.background), rescale.p=TRUE)
      r["chisq"] = a$statistic
      r["p"] = a$p.value
    }, error = function(e) r["p"] = 1,
    warning = function(e) r["p"] = 1)
  }

  r
}

# Looks for enrichment of a motif upstream of a group of genes.
# Args:
#   motif.count - motif counts at various cutoffs
#   cluster - genes in group, as a logical vector indexed by gene name
# Returns: array of results
motif.enriched.in.group = function(motif.count, cl) {

  g = intersect(names(cl), dimnames(motif.count)[[1]])

  # counts of motifs
  motifs.1 = motif.count[g,,,]
  motifs.cl = apply(motifs.1[cl,,,], c(2,3,4), sum)
  motifs.bg = apply(motifs.1[!cl,,,], c(2,3,4), sum)

  # amount of sequence being considered
  upstream.seq.1 = upstream.region.size[g,,]
  upstream.cl = apply(upstream.seq.1[cl,,], c(2,3), sum)
  upstream.bg = apply(upstream.seq.1[!cl,,], c(2,3), sum)

  # array of results
  r = array.from.dimnames(list(
    upstream.dist.kb = dimnames(motif.count)[[2]],
    conservation = dimnames(motif.count)[[3]],
    motif.score = dimnames(motif.count)[[4]],
    stat = stat.names))

  for(ud in dimnames(motifs.1)$upstream.dist)
    for(cn in dimnames(motifs.1)$conservation)
      for(sc in dimnames(motifs.1)$score) {
        r[ ud, cn, sc, ] =
          compute.enrichment(motifs.cl[ud, cn, sc],
            motifs.bg[ud, cn, sc],
            upstream.cl[ud, cn],
            upstream.bg[ud, cn])
      }

  r
}

# Looks for motifs enriched upstream of various groups of genes.
# Args:
#   motif.count.dir - directory containing motif counts (this will
#     search all motifs in that directory)
#   cluster - the clusters in which to search for enrichment,
#     as a list of boolean vectors, each indexed by gene name
#   motifs - (optional) which motifs to check for enrichment
# Returns: array with indices "cluster", "motif", and "cutoff",
#   plus "stat", consisting of "motifs.cluster", "motifs.background",
#   "bp.cluster", "bp.background", "chisq", "p", and "p.corr".
motifs.enriched.in.groups = function(motif.count.dir, clusters, motifs=NULL) {
  motifs.in.dir = sub(".Rdata", "", list.files(motif.count.dir))

  if (is.null(motifs)) {
    motifs = motifs.in.dir
  }
  else {
    motifs = intersect(motifs, motifs.in.dir)
  }

  # XXX hack to get dimension names from one of the motif counts
  motif.count = NULL
  load(paste0(motif.count.dir, "/", motifs[1], ".Rdata"))

  # array of results
  r = array.from.dimnames(list(motif = motifs,
    group = as.character(names(clusters)),
    upstream.dist.kb = dimnames(motif.count)[[2]],
    conservation = dimnames(motif.count)[[3]],
    motif.score = dimnames(motif.count)[[4]],
    stat = stat.names))

  # loop through the motifs
  for(motif in motifs) {

    # read in motif counts
    motif.count = NULL
    load(paste0(motif.count.dir, "/", motif, ".Rdata"))

    # loop through the clusters
    for(cl in names(clusters)) {

      write.status(paste(motif, cl))
      r[motif, cl,,,, ] =
        motif.enriched.in.group(motif.count, clusters[[cl]])
      gc()
    }
  }

  r
}

# Does a large number of chi-square tests.
# All its args can be vectors of the same length; you may want to
# reshape from array to vector, and back, to use these.
# Args:
#   motifs.cl, motifs.bg - # of motifs in cluster and background
#   bp.cl, bp.bg - # of base pairs in cluster and background
# Returns: matrix with one row per entry of each of these, and
#   one column for each of the statistics
chisq.test.many = function(motifs.cl, motifs.bg, bp.cl, bp.bg) {
  n = length(motifs.cl)

  # array of results
  r = matrix(nrow=n, ncol=8)
  colnames(r) = stat.names
  r[ , "motifs.cluster" ] = motifs.cl
  r[ , "motifs.background" ] = motifs.bg
  r[ , "bp.cluster" ] = bp.cl
  r[ , "bp.background" ] = bp.bg

  motifs = motifs.cl + motifs.bg
  bp = bp.cl + bp.bg

  # expected number of motifs in cluster, and background
  e = motifs * (bp.cl / bp)
  e.bg = motifs * (bp.bg / bp)

  r[ , "enrich" ] = motifs.cl / e
  r[ , "chisq" ] = ((motifs.cl - e) ^ 2) / e +
    ((motifs.bg - e.bg) ^ 2) / e.bg
  r[ , "p" ] = pchisq(r[ , "chisq" ], 1, lower.tail=FALSE)
  
  r
}

# Faster version of the above (or an attempt at this.)
# Args:
#   cl - the clustering, as a vector of integers (or strings),
#     indexed by gene
#   motif.count - array of motif counts
#   upstream.bp - array of upstream counts
# Returns: array of results
motifs.enriched.in.group.2 = function(cl, motif.count, upstream.bp) {
  cl.names = names(cl)
  cl = as.character(cl)
  names(cl) = cl.names

  clusters = sort(unique(cl))

  motif.count = motif.count[names(cl),,,]
  upstream.bp = upstream.bp[names(cl),,]

  # count motifs by cluster (foreground and background)
  motifs.cluster =
    aperm(apply(motif.count, c(2,3,4), function(x) by(x, cl, sum)), c(2,3,4,1))
  motifs.total = apply(motifs.cluster, c(1,2,3), sum)
  # aperm() was used so that this works
  motifs.background = as.vector(motifs.total) - motifs.cluster

  # add up amount of upstream sequence (foreground and background)
  bp.cluster = aperm(apply(upstream.bp, c(2,3),
    function(x) by(x, cl, sum)), c(2,3,1))
  bp.total = apply(bp.cluster, c(1,2), sum)
  bp.background = as.vector(bp.total) - bp.cluster

  # compute chi-squared scores (again, using aperm() gymnastics)
  r1 = chisq.test.many(as.vector(aperm(motifs.cluster, c(4,1,2,3))),
    as.vector(aperm(motifs.background, c(4,1,2,3))),
    rep(as.vector(aperm(bp.cluster, c(3,1,2))), 3),
    rep(as.vector(aperm(bp.background, c(3,1,2))), 3))

  list(motifs.cluster = motifs.cluster,
    motifs.background = motifs.background,
    bp.cluster = bp.cluster, bp.background = bp.background, r1 = r1)
}

