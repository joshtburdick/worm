# Looks for motifs and ChIP signals enriched upstream
# of genes in clusters (or other groups.)
# Optimizes cutoffs to maximize number of things
# significantly enriched.
# Note that I'm calling ChIP signals "motifs" as well, which
# is somewhat inaccurate.

# FIXME
# - alter this to avoid "chi-squared may be inaccurate" warnings?
# - write out .tsv files directly, instead of .Rdata files?

source("git/utils.r")
source("git/sort_paper/enrichment/groupToList.r")

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

# XXX for testing
# known.motifs.small = "RFX2_DBD"

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



# Does a chi-squared test for a clustering.
# Args:
#   gene.in.cluster - the clustering to use, as a logical vector indexed by gene name
#   motif.count - the number of motif occurences (as a numeric vector,
#     indexed by gene name)
#   upstream.size - the size of the upstream regions for each gene
#     (as a numeric vector indexed by gene)
# Returns: a vector of statistics
compute.enrichment.1 = function(gene.in.cluster, motif.count, upstream.size) {

  # XXX note some hackery to deal with different gene names
#  cl1 = intersect(names(gene.in.cluster), names(motif.count))
#  mc1 = motif.count[cl1]
#  motifs.cluster = sum(mc1[ cl1 ])
#  motifs.background = sum(mc1[ !cl1 ])

  motifs.cluster = sum(motif.count[ gene.in.cluster ], na.rm=TRUE)
  motifs.background = sum(motif.count[ ! gene.in.cluster ], na.rm=TRUE)

  bpc = upstream.size[ names(gene.in.cluster) ]
  bp.cluster = sum( bpc[ gene.in.cluster ], na.rm=TRUE)
  bp.background = sum( bpc[ ! gene.in.cluster ], na.rm=TRUE)

  # compute counts
  r = c(motifs.cluster = motifs.cluster,
    motifs.background = motifs.background,
    bp.cluster = bp.cluster,
    bp.background = bp.background,
    enrichment =
    (motifs.cluster / motifs.background) / (bp.cluster / bp.background),
    chisq = NA, p = 1, p.corr = NA)

# browser()   # XXX for testing

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

# Looks for enrichments with various cutoffs.
# Args:
#   motifs - list of "motifs" (which are either TFs or ChIP signals)
#   get.counts - function which gets counts for a motif with a given name
#   clusters - the clusters to use, as a list of logical vectors
#   upstream.dist.kb.cutoff - the cutoffs to use for upstream distance
#   conservation.cutoff -                      " conservation
#   score.cutoff -                             " for the "score" column
# Returns: array with indices as above, plus
#     "motifs.cluster", "motifs.background", "bp.cluster", "bp.background",
#     "chisq", "p", and "p.corr".
compute.enrichment.diff.cutoffs = function(motifs, get.counts, clusters,
    upstream.dist.kb.cutoff, conservation.cutoff, score.cutoff) {

  # array of results
  stat.names = c("motifs.cluster", "motifs.background",
    "bp.cluster", "bp.background", "enrich", "chisq", "p", "p.corr")
  r = array.from.dimnames(list(motif = motifs,
    group = as.character(names(clusters)),
    upstream.dist.kb = upstream.dist.kb.cutoff,
    conservation = conservation.cutoff,
    motif.score = score.cutoff,
    stat = stat.names))

  # cache of upstream conservation at different cutoffs,
  # indexed first by upstream distance, then conservation
  upstream.cons.at.cutoff = NULL
  for(upstream.dist.kb in 1:5) {
    for(conservation in c(0, 0.5, 0.7, 0.9)) {
      bp1 = upstream.cons.dist[[upstream.dist.kb]]
      upstream.bp =
        apply(bp1[ , colnames(bp1) >= conservation ], 1, sum)
      upstream.bp = upstream.bp[ names(clusters[[1]]) ]
      upstream.cons.at.cutoff[[as.character(upstream.dist.kb)]][[as.character(conservation)]] = upstream.bp
    }
  }

  # loop through the motifs
  for (motif in motifs) {
    m = get.motif.counts(motif)

    # loop through various cutoffs
    for (upstream.dist.kb in upstream.dist.kb.cutoff)
      for (conservation in conservation.cutoff) {

        # amount of upstream sequence present at that cutoff
        # FIXME cache this
if (FALSE) {
        bp1 = upstream.cons.dist[[upstream.dist.kb]]
        upstream.bp = apply(bp1[ , colnames(bp1) >= conservation ], 1, sum)
        upstream.bp = upstream.bp[ names(clusters[[1]]) ]
}
        upstream.bp = upstream.cons.at.cutoff[[as.character(upstream.dist.kb)]][[as.character(conservation)]]

        # loop through motif score cutoffs, and count motifs at that cutoff
        for (motif.score in score.cutoff) {
          write.status(paste(motif, upstream.dist.kb,
            conservation, motif.score))

          m1 = m[ m$upstream.dist >= -1000 * upstream.dist.kb &
            m$motif.cons >= conservation &
            m$motif.score >= motif.score , ]
          mc1 = c(table(m1$gene))
          motif.counts = rep(0, length(names(clusters[[1]])))
          motif.counts[ names(mc1) ] = mc1
          gc()

          for(cl in as.character(names(clusters))) {
            r[motif,cl,as.character(upstream.dist.kb),
              as.character(conservation),
              as.character(motif.score),] =
                compute.enrichment.1(clusters[[cl]], motif.counts, upstream.bp)
          }
        }
      }
  }

  r
}

# Runs this on clustering for the motifs.
find.cluster.enrichment.motifs = function() {
  output.dir = "git/sort_paper/tf/motif/cluster.enrichment/"
  system(paste("mkdir -p", output.dir))

  for (f in list.files(clustering.dir)) {
    cat(paste(f, date), "\n")

    # clustering to use
#    clustering1 = read.tsv(paste0(clustering.dir, "/", f, "/clusters.tsv"))
#    clustering = clustering1[,2]
#    names(clustering) = rownames(clustering1)

    clusters = cluster.to.gene.list(paste0(clustering.dir, "/", f, "/clusters.tsv"))

    enrich = compute.enrichment.diff.cutoffs(known.motifs.small, get.motif.counts,
      clusters, c(1,2,3), c(0, 0.5, 0.7, 0.9), c(30,35,40))
    enrich[,,,,,"p"][ is.na(enrich[,,,,,"p"]) ] = 1

    enrich[,,,,,"p.corr"] = array(p.adjust(as.vector(enrich[,,,,,"p"]), method="fdr"),
      dim=dim(enrich[,,,,,"p.corr"]), dimnames=dimnames(enrich[,,,,,"p.corr"]))

    save(enrich, file=paste0(output.dir, "/", f, ".Rdata"))
  }
}

find.cluster.enrichment.motifs()

