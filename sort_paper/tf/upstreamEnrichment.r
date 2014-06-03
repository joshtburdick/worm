# Looks for motifs and ChIP signals enriched upstream
# of genes in clusters (or other groups.)
# Optimizes cutoffs to maximize number of things
# significantly enriched.
# Note that I'm calling ChIP signals "motifs" as well, which
# is somewhat inaccurate.

# FIXME: alter this to avoid "chi-squared may be inaccurate"
# warnings?

source("git/utils.r")

clustering.dir = "git/cluster/hierarchical/"

# running this on the filtered (non-redundant) motifs
motif.gene.dir = "git/cluster/motif/distAndConservation/5kb/"
known.motifs = {
  r = list.files(motif.gene.dir)
  sub("_upstreamMotifCons.tsv.gz", "", r)
}
motif.filter = read.tsv("git/tf/motif/motifFilter.tsv")
known.motifs.small =
  intersect(known.motifs, motif.filter$canonical.name)

# Gets the counts of one motif upstream of all genes.
get.motif.counts = function(m) {

  r = read.table(paste(motif.gene.dir, m, "_upstreamMotifCons.tsv.gz", sep=""),
    as.is=TRUE)
  colnames(r) =
    c("region.chr", "region.a", "region.b",
    "gene", "score", "strand",
    "motif", "motif.a", "motif.b", "motif.id", "motif.score",
    "motif.strand", "motif.cons")
  r$upstream.dist = ifelse(r$strand=="+",
    r$motif.b - r$region.b, r$region.a - r$motif.a)
  r
}

# number of upstream bp with different levels of conservation
upstream.cons.dist = list()
for(i in 1:5)
  upstream.cons.dist[[i]] =
    read.tsv(paste0("git/tf/motif/conservation/cons_hist_WS220_",
      i, "kb_upstream.tsv.gz"))

# Does a chi-squared test for each of a set of clusters.
# Args:
#   motif.name - name of the motif
#   clustering - the clustering to use, as a numeric vector
#     indexed by gene name
#   m - the number of motif occurences (as a numeric vector,
#     indexed by cluster name)
#   upstream.size - the size of the upstream regions for each
#     gene (as a vector indexed by gene)
# Returns: vector of p-values
compute.enrichment.1 = function(motif.name, clustering, m, upstream.size) {
  clustering = clustering[ !is.na(clustering) ]

  clusters = as.character(sort(unique(clustering)))
  p = rep(NA, length(clusters))
  names(p) = clusters

  # make sure gene names match
  g = intersect(names(clustering), names(upstream.size))
  clustering = clustering[ g ]
  upstream.size = upstream.size[ g ]

  # total region size for each cluster
  bp.per.cluster = c(by(upstream.size, clustering, sum))

  # count of motifs for each cluster
  motifs.per.cluster = rep(0, length(clusters))
  names(motifs.per.cluster) = clusters
  motifs.per.cluster[names(m)] = m

  motifs.total = sum(motifs.per.cluster)
  bp.total = sum(bp.per.cluster)

  # array to hold results
  stat.names = c("motifs.cluster", "motifs.background",
    "bp.cluster", "bp.background", "chisq", "p", "p.corr")
  r = array.from.dimnames(list(group=clusters, stat=stat.names))

  for(cl in as.character(clusters)) {

    # compute counts
    r[cl,"motifs.cluster"] = motifs.per.cluster[cl]
    r[cl,"motifs.background"] = motifs.total - r[cl,"motifs.cluster"]
    r[cl,"bp.cluster"] = bp.per.cluster[cl]
    r[cl,"bp.background"] = bp.total - r[cl,"bp.cluster"]

    # do the chi-square test; on error, just set the p-value to 1
    tryCatch({
      a = chisq.test(c(r[cl,"motifs.cluster"], r[cl,"motifs.background"]),
        p = c(r[cl,"bp.cluster"], r[cl,"bp.background"]), rescale.p=TRUE)
        r[cl,"chisq"] = a$statistic
        r[cl,"p"] = a$p.value
    }, error = function(e) r[cl,"p"] = 1,
      warning = function(e) r[cl,"p"] = 1)

if (FALSE) {
    if (abs(r[cl,"motifs.cluster"] + r[cl,"motifs.background"]) > 0) {

      # ??? should this be one-sided?

      a = chisq.test(c(r[cl,"motifs.cluster"], r[cl,"motifs.background"]),
        p = c(r[cl,"bp.cluster"], r[cl,"bp.background"]), rescale.p=TRUE)
      # only include enrichments
      if (a$observed[1] > a$expected[1]) {
        r[cl,"chisq"] = a$statistic
        r[cl,"p"] = a$p.value
      }
      else {
        r[cl,"p"] = 1
      }
    }
    else {
      r[cl,"p"] = 1
    }
}

  }

  r
}

# Faster version of the above (hopefully.)
# Args:
#   motifs - list of "motifs" (which are either TFs or ChIP signals)
#   get.counts - function which gets counts for a motif with a given name
#   clustering - the clustering to use
#   upstream.dist.kb.cutoff - the cutoffs to use for upstream distance
#   conservation.cutoff -                      " conservation
#   score.cutoff -                             " for the "score" column
# Returns: array with indices as above, plus
#     "motifs.cluster", "motifs.background", "bp.cluster", "bp.background",
#     "chisq", "p", and "p.corr".
compute.enrichment.diff.cutoffs = function(motifs, get.counts, clustering,
    upstream.dist.kb.cutoff, conservation.cutoff, score.cutoff) {
  clustering = clustering[ !is.na(clustering) ]
  clusters = as.character(sort(unique(clustering)))

  # array of results
  stat.names = c("motifs.cluster", "motifs.background",
    "bp.cluster", "bp.background", "chisq", "p", "p.corr")
  r = array.from.dimnames(list(motif = motifs,
    group = clusters,
    upstream.dist.kb = upstream.dist.kb.cutoff,
    conservation = conservation.cutoff,
    motif.score = score.cutoff,
    stat = stat.names))

  # loop through the motifs
  for (motif in motifs) {
    m = get.motif.counts(motif)

    # try various cutoffs
    for (upstream.dist.kb in upstream.dist.kb.cutoff)
      for (conservation in conservation.cutoff) {

        # amount of upstream sequence present at that cutoff
        bp1 = upstream.cons.dist[[upstream.dist.kb]]
        upstream.bp = apply(bp1[ , colnames(bp1) >= conservation ], 1, sum)

        for (motif.score in score.cutoff) {
          write.status(paste(motif, upstream.dist.kb,
            conservation, motif.score))

          m1 = m[ m$upstream.dist >= -1000 * upstream.dist.kb &
            m$motif.cons >= conservation &
            m$motif.score >= motif.score , ]
          gc()
          m.counts = c(table(clustering[m1$gene]))

          en = compute.enrichment.1(motif, clustering, m.counts, upstream.bp)

          r[motif,,as.character(upstream.dist.kb),as.character(conservation),as.character(motif.score),] =
            en[clusters,]
        }
      }
  }

  r
}

# Runs this on clustering for the motifs.
find.cluster.enrichment.motifs = function() {
  output.dir = "git/sort_paper/tf/motif/"
  system(paste("mkdir -p", output.dir))

  for (f in list.files(clustering.dir)) {
    cat(f, "\n")

    # clustering to use
    clustering1 = read.tsv(paste0(clustering.dir, "/", f, "/clusters.tsv"))
    clustering = clustering1[,2]
    names(clustering) = rownames(clustering1)

    enrich = compute.enrichment.diff.cutoffs(known.motifs.small, get.motif.counts,
      clustering, c(1,2,3), c(0, 0.5, 0.7, 0.9), c(30,35,40))
    enrich[,,,,,"p"][ is.na(enrich[,,,,,"p"]) ] = 1

    enrich[,,,,,"p.corr"] = array(p.adjust(as.vector(enrich[,,,,,"p"]), method="fdr"),
      dim=dim(enrich[,,,,,"p.corr"]), dimnames=dimnames(enrich[,,,,,"p.corr"]))

    save(enrich, file=paste0(output.dir, "/", f, ".Rdata"))
  }
}

# find.cluster.enrichment.motifs()

