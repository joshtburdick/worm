# Optimizes cutoffs to maximize number of motifs
# significantly enriched.

# FIXME: alter this to avoid "chi-squared may be inaccurate"
# warnings?

source("git/utils.r")

clustering.dir = "git/cluster/hierarchical/"

motif.gene.dir = "git/cluster/motif/distAndConservation/5kb/"
output.dir = "git/cluster/motif/enrichOptimize/cutoff.optimize"
system(paste("mkdir -p", output.dir))

# running this on the filtered (non-redundant) motifs
known.motifs = {
  r = list.files(motif.gene.dir)
  sub("_upstreamMotifCons.tsv.gz", "", r)
}
motif.filter = read.tsv("git/tf/motif/motifFilter.tsv")
known.motifs.small =
  intersect(known.motifs, motif.filter$canonical.name)

# number of upstream bp with different levels of conservation
upstream.cons.dist = list()
for(i in 1:5)
  upstream.cons.dist[[i]] =
    read.tsv(paste0("git/tf/motif/conservation/cons_hist_WS220_",
      i, "kb_upstream.tsv.gz"))

# clustering to use
clustering1 = read.tsv(
  "git/cluster/hierarchical/hier.50.clusters/clusters.tsv")
clustering = clustering1[,2]
names(clustering) = rownames(clustering1)

# Gets the motif counts for one gene.
get.motif.counts = function(m) {

  r = read.table(paste(motif.gene.dir, m, "_upstreamMotifCons.tsv.gz", sep=""),
    as.is=TRUE)
  colnames(r) =
    c("region.chr", "region.a", "region.b",
    "gene", "score", "strand",
    "motif", "motif.a", "motif.b", "motif.id", "motif.score",
    "motif.strand", "motif.cons")
  r$upstream.dist = ifelse(r$strand=="+",
    r$motif.a - r$region.b, r$region.a - r$motif.a)
  r
}

# Computes chi-squared enrichment stats for a given set of
# clusters, and a set of motifs.
# Args:
#   motif.name - name of the motif
#   clustering - the clustering to use, as a numeric vector
#     indexed by gene name
#   m - the number of motif occurences (as a numeric vector,
#     indexed by cluster name)
#   upstream.size - the size of the upstream regions for each
#     gene (as a vector indexed by gene)
compute.enrichment.1 = function(motif.name, clustering, m, upstream.size) {
  clustering = clustering[ !is.na(clustering) ]

  clusters = sort(unique(clustering))
  r = NULL

# print(clustering[1:10])
# print(upstream.size[1:10])

  # make sure gene names match
  g = intersect(names(clustering), names(upstream.size))
# print(g[1:20])
  clustering = clustering[ g ]
  upstream.size = upstream.size[ g ]

  # summaries for each cluster
# print(upstream.size[1:20])
# print(clustering[1:20])

  bp.per.cluster = c(by(upstream.size, clustering, sum))

  # count of motifs for each cluster
  motifs.per.cluster = rep(0, length(clusters))
  names(motifs.per.cluster) = clusters
  motifs.per.cluster[names(m)] = m

  motifs.total = sum(motifs.per.cluster)
  bp.total = sum(bp.per.cluster)

  for(cl in as.character(clusters)) {

#    write.status(cl)

    # compute counts
    motifs.cluster = motifs.per.cluster[cl]
    motifs.background = motifs.total - motifs.cluster
    bp.cluster = bp.per.cluster[cl]
    bp.background = bp.total - bp.cluster

    # ??? should this be one-sided?
    a = chisq.test(c(motifs.cluster, motifs.background),
      p = c(bp.cluster, bp.background), rescale.p=TRUE)
    r1 = data.frame(motif = motif.name, cluster = cl,
      motifs.cluster = motifs.cluster, motifs.background = motifs.background,
      bp.cluster = bp.cluster, bp.background = bp.background,
      enrichment = (motifs.cluster / bp.cluster) /
        (motifs.background / bp.background),
      X.squared = a$statistic, p = a$p.value, stringsAsFactors=FALSE)

    r = rbind(r, r1)
  }

  r
}

# Like the above, but just returns a vector of unadjusted p-values.
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

  # count of motifs for each cluster
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

# Computes enrichment with several different cutoffs
# for various things.
compute.enrichment.diff.cutoffs = function(motif, clustering) {

  m = get.motif.counts(motif)

  # try various cutoffs
  for (upstream.dist.kb in c(1:3))
    for (conservation in c(0.5, 0.7, 0.9)) {

      # amount of upstream sequence present at that cutoff
      bp1 = upstream.cons.dist[[upstream.dist.kb]]
      upstream.bp = apply(bp1[ , colnames(bp1) >= conservation ], 1, sum)

      for (motif.score in c(30, 35, 40)) {

        write.status(paste(motif, upstream.dist.kb,
          conservation, motif.score))

        m1 = m[ m$upstream.dist >= -1000 * upstream.dist.kb &
          m$motif.cons >= conservation &
          m$motif.score >= motif.score , ]
        gc()
        m.counts = c(table(clustering[m1$gene]))

        r1 = compute.enrichment.1(motif, clustering, m.counts, upstream.bp)
        r = rbind(r, data.frame(r1, stringsAsFactors=FALSE))
      }
    }
  r
}

# Faster version of the above (hopefully.)
# Args:
#   motifs - list of motifs
#   clustering - the clustering to use
# Returns: a large array of unadjusted p-values
compute.enrichment.diff.cutoffs.faster = function(motifs, clustering) {
  clustering = clustering[ !is.na(clustering) ]
  clusters = as.character(sort(unique(clustering)))

  # array of p-values
  r = array(dim = list(length(motifs), length(clusters), 3, 4, 3),
    dimnames = list(motif = motifs,
      group = clusters,
      upstream.dist.kb = c(1, 2, 3),
      conservation = c(0, 0.5, 0.7, 0.9),
      motif.score = c(30, 35, 40)))

  # loop through the motifs
  for (motif in motifs) {
    m = get.motif.counts(motif)

    # try various cutoffs
    for (upstream.dist.kb in c(1:3))
      for (conservation in c(0, 0.5, 0.7, 0.9)) {

        # amount of upstream sequence present at that cutoff
        bp1 = upstream.cons.dist[[upstream.dist.kb]]
        upstream.bp = apply(bp1[ , colnames(bp1) >= conservation ], 1, sum)

        for (motif.score in c(30, 35, 40)) {
          write.status(paste(motif, upstream.dist.kb,
            conservation, motif.score))

          m1 = m[ m$upstream.dist >= -1000 * upstream.dist.kb &
            m$motif.cons >= conservation &
            m$motif.score >= motif.score , ]
          gc()
          m.counts = c(table(clustering[m1$gene]))

          p = compute.enrichment.p.only(motif, clustering, m.counts, upstream.bp)
          r[motif,,as.character(upstream.dist.kb),as.character(conservation),as.character(motif.score)] =
            p[clusters]
        }
      }
  }

  r
}

if (FALSE) {
for (f in list.files(clustering.dir)[1]) {
  cat(f, "\n")

  # clustering to use
  clustering1 = read.tsv(paste0(clustering.dir, "/", f, "/clusters.tsv"))
  clustering = clustering1[,2]
  names(clustering) = rownames(clustering1)

  r = NULL
  for(m in known.motifs.small) {
    r = rbind(r, compute.enrichment.diff.cutoffs(m, clustering))
  }

  r$p.corr = p.adjust(r$p, method="fdr")
  r = r[ r$p.corr <= 0.05 & r$X.squared > 0 , ]

  write.tsv(r, paste0(output.dir, "/", f, ".tsv"))
}
}

if (TRUE) {
for (f in list.files(clustering.dir)) {
  cat(f, "\n")

  # clustering to use
  clustering1 = read.tsv(paste0(clustering.dir, "/", f, "/clusters.tsv"))
  clustering = clustering1[,2]
  names(clustering) = rownames(clustering1)

  p = compute.enrichment.diff.cutoffs.faster(
    known.motifs.small, clustering)
  p.corr = array(p.adjust(as.vector(p), method="fdr"),
    dim=dim(p), dimnames=dimnames(p))

  save(p, p.corr, file=paste0(output.dir, "/", f, ".Rdata"))
}
}

