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

