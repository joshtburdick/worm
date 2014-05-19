# Optimizes cutoffs to maximize number of ChIP peaks
# significantly enriched (based on similar motif thing.)

# FIXME: alter this to avoid "chi-squared may be inaccurate"
# warnings?

source("git/utils.r")

clustering.dir = "git/cluster/hierarchical/"

chip.gene.dir = "git/cluster/chip/distAndConservation/"
output.dir = "git/cluster/chip/enrichOptimize/cutoff.optimize"

system(paste("mkdir -p", output.dir))

# names of experiments
chip.names = list.files("git/cluster/chip/distAndConservation")
chip.names = grep("rep1", chip.names, value=TRUE)
chip.names = sub("_upstreamChipCons.tsv.gz", "", chip.names)

# number of upstream bp with different levels of conservation
upstream.cons.dist = list()
for(i in 1:5)
  upstream.cons.dist[[i]] =
    read.tsv(paste0("git/tf/motif/conservation/cons_hist_WS220_",
      i, "kb_upstream.tsv.gz"))

# clustering to use (deprecated)
clustering1 = read.tsv(
  "git/cluster/hierarchical/hier.50.clusters/clusters.tsv")
clustering = clustering1[,2]
names(clustering) = rownames(clustering1)

# Gets the motif counts for one gene.
get.chip.counts = function(m) {

  r = read.table(paste(chip.gene.dir, m, "_upstreamChipCons.tsv.gz", sep=""),
    as.is=TRUE)
  colnames(r) =
    c("region.chr", "region.a", "region.b",
    "gene", "score", "strand",
    "chip.chr", "chip.a", "chip.b", "chip.id", "chip.score",
    "chip.strand", "chip.cons")
  r$upstream.dist = ifelse(r$strand=="+",
    r$chip.b - r$region.b, r$region.a - r$chip.a)
  r
}

# Like the above, but just returns a vector of unadjusted p-values
# (and is faster.)
# Args:
#   motif.name - name of the motif
#   clustering - the clustering to use, as a numeric vector
#     indexed by gene name
#   m - the number of motif occurences (as a numeric vector,
#     indexed by cluster name)
#   upstream.size - the size of the upstream regions for each
#     gene (as a vector indexed by gene)
# Returns: vector of p-values
compute.enrichment.p.only = function(motif.name, clustering, m, upstream.size) {
  clustering = clustering[ !is.na(clustering) ]

  clusters = sort(unique(clustering))
  p = rep(NA, length(clusters))
  names(p) = clusters

  # make sure gene names match
  g = intersect(names(clustering), names(upstream.size))
  clustering = clustering[ g ]
  upstream.size = upstream.size[ g ]

  # total region size for each cluster
  bp.per.cluster = c(by(upstream.size, clustering, sum))

  # count of motifs (er, ChIP peaks) for each cluster
  motifs.per.cluster = rep(0, length(clusters))
  names(motifs.per.cluster) = clusters
  motifs.per.cluster[names(m)] = m

  motifs.total = sum(motifs.per.cluster)
  bp.total = sum(bp.per.cluster)

  for(cl in as.character(clusters)) {

    # compute counts
    motifs.cluster = motifs.per.cluster[cl]
    motifs.background = motifs.total - motifs.cluster
    bp.cluster = bp.per.cluster[cl]
    bp.background = bp.total - bp.cluster

    if (abs(motifs.cluster + motifs.background) > 0) {
      # ??? should this be one-sided?
      a = chisq.test(c(motifs.cluster, motifs.background),
        p = c(bp.cluster, bp.background), rescale.p=TRUE)
      # only include enrichments
      if (a$observed[1] > a$expected[1]) {
        p[cl] = a$p.value
      }
      else {
        p[cl] = 1
      }
    }
    else {
      p[cl] = 1
    }
  }

  p
}

# Computes enrichments for many possible values of cutoffs.
# Args:
#   chip.names - name of the chip experiments
#   clustering - the clustering to use
# Returns: a large array of unadjusted p-values
compute.enrichment.diff.cutoffs.faster = function(chip.names, clustering) {
  clustering = clustering[ !is.na(clustering) ]
  clusters = as.character(sort(unique(clustering)))

  # array of p-values
  r = array(dim = list(length(chip.names), length(clusters), 3, 4, 4),
    dimnames = list(chip = chip.names,
      group = clusters,
      upstream.dist.kb = c(1, 2, 3),
      conservation = c(0, 0.5, 0.7, 0.9),
      chip.score.quantile = c(0, 0.25, 0.5, 0.75)))

  # loop through the motifs
  for (chip in chip.names) {
    m = get.chip.counts(chip)

    # try various cutoffs
    for (upstream.dist.kb in c(1:3))
      for (conservation in c(0, 0.5, 0.7, 0.9)) {

        # amount of upstream sequence present at that cutoff
        bp1 = upstream.cons.dist[[upstream.dist.kb]]
        upstream.bp = apply(bp1[ , colnames(bp1) >= conservation ], 1, sum)

        for (chip.score.quantile in c(0, 0.25, 0.5, 0.75)) {
          write.status(paste(chip, upstream.dist.kb,
            conservation, chip.score.quantile))

          # compute cutoff of chip score to use
          chip.score.cutoff = quantile(m$chip.score, prob=chip.score.quantile)

          m1 = m[ m$upstream.dist >= -1000 * upstream.dist.kb &
            m$chip.cons >= conservation &
            m$chip.score >= chip.score.cutoff , ]
          gc()
          m.counts = c(table(clustering[m1$gene]))

          p = compute.enrichment.p.only(chip, clustering, m.counts, upstream.bp)
          r[chip,,as.character(upstream.dist.kb),as.character(conservation),as.character(chip.score.quantile)] =
            p[clusters]
        }
      }
  }

  r
}

if (TRUE) {
for (f in list.files(clustering.dir)) {
  cat(f, "\n")

  # clustering to use
  clustering1 = read.tsv(paste0(clustering.dir, "/", f, "/clusters.tsv"))
  clustering = clustering1[,2]
  names(clustering) = rownames(clustering1)

  p = compute.enrichment.diff.cutoffs.faster(
    chip.names, clustering)
  p.corr = array(p.adjust(as.vector(p), method="fdr"),
    dim=dim(p), dimnames=dimnames(p))

  save(p, p.corr, file=paste0(output.dir, "/", f, ".Rdata"))
}
}

