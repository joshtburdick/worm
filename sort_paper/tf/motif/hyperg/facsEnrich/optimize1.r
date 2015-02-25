# Attempt at comparing the enrichment of motifs in
# FACS-sorted fractions, with different background models.

# XXX
options(stringsAsFactors = FALSE)

library("reshape2")

source("git/sort_paper/tf/motif/hyperg/hypergTestGroups.r")
source("git/utils.r")

# names of shuffled motifs
shuffled.motifs = sub(".Rdata", "",
  list.files("git/tf/motif/count/upstreamMotifCount/jolma2013_shuffled/"))

motifs = intersect(orig.motif.list, shuffled.motifs)

# for testing
foxb.motifs = c("FOXB1_DBD", "FOXB1_DBD_2", "FOXB1_DBD_3", "FOXB1_full")

# Max. expression in reads per million
rpm = read.tsv("git/cluster/readsPerMillion.tsv")
rpm = rpm[ , !grepl("^(HS|RNAi)", colnames(rpm)) ]
max.expr = apply(rpm, 1, max)

# "Bucketizes" expression at some level.
bucketize = function(x) floor(log2(1+x))

# "Buckets" of genes expressed at a similar level.
expr.by.bucket = by(names(max.expr), bucketize(max.expr),
  function(x) c(as.character(x)))

# Creates a set of genes with matched expression levels.
# Args:
#   g - a list of genes
# Returns: a list of genes with similar expression levels.
#   (The list may be shorter, if some of the buckets don't
#   contain enough genes.)
matched.expr.genes = function(g) {
  r = NULL

  buckets = bucketize(max.expr[g])
  for(b in unique(buckets)) {
    g1 = names(buckets)[ buckets == b ]
    a = setdiff(expr.by.bucket[[as.character(b)]], g1)
    if (length(a) >= length(g1))
      r = c(r, sample(a, length(g1)))
    else
      r = c(r, a)
  }

  r
}

# Quick test of matched.expr.genes().
matched.expr.genes.test1 = function(n) {
  par(mfrow=c(1,2))

  g = sample(names(max.expr), n)
  hist(log2(1+max.expr[g]), breaks=100, col="grey")

  g1 = matched.expr.genes(g)
  hist(log2(1+max.expr[g1]), breaks=100, col="grey")
}

# Get enrichments for motifs, and shuffled motifs.
# Args:
#   s - a list of gene sets (as boolean vectors)
# Returns: enrichment for original, and shuffled motifs.
enrich.test = function(s) {

  # only keep cases in which at least some genes are in
  # both the "foreground" and "background" sets
  s = s[ sapply(s, var) > 0 ]

  r1 = enrich.test.gene.sets.many.motifs(
    "git/tf/motif/count/upstreamMotifCount/5kb/",
    s, motifs)    # was foxb.motifs
  r = melt(r1[,,"p",3,3,3])

  r1.shuffled = enrich.test.gene.sets.many.motifs(
    "git/tf/motif/count/upstreamMotifCount/jolma2013_shuffled/",
    s, motifs)   # was foxb.motifs)
  r.shuffled = melt(r1.shuffled[,,"p",3,3,3])

  colnames(r)[[2]] = "sort_marker"
  colnames(r)[[3]] = "p"
  # XXX this assumes things are in the same order in r and r.shuffled
  cbind(r, shuffled.p = r.shuffled[,3],
    stringsAsFactors=FALSE)
}

# Sets of genes enriched, versus genes with no change
# (this should be the same as previous.)
# Args:
#   r - a data table
#   cutoff - the cutoff to use
# Returns: a list of boolean vectors, one per column of r.
enriched.vs.no.change = function(r, cutoff) {
  a = NULL
  for(s in colnames(r)) {
    a[[ s ]] = (r[,s] >= cutoff)[ r[,s] > -cutoff ]
  }
  a
}

# Sets of genes enriched, versus depleted.
# Args:
#   r - a data table
#   cutoff - the cutoff to use
# Returns: a list of boolean vectors, one per column of r.
enriched.vs.depleted = function(r, cutoff) {
  a = NULL
  for(s in colnames(r)) {
    a[[ s ]] =
      (r[,s] >= cutoff)[ abs(r[,s]) > cutoff ]
  }
  a
}

# Utility function.
tf.vector = function(true, false) {
  a = c(true, false)
  r = a %in% true
  names(r) = a
  r
}

# Sets of genes enriched (or depleted), versus a random set.
# Args:
#   r - a data table
#   cutoff - the cutoff to use
# Returns: a list of boolean vectors, one per column of r.
enriched.versus.random = function(r, cutoff) {
  a = NULL
  for(s in colnames(r)) {
    g = rownames(r)[ r[,s] >= cutoff ]
    g.random = matched.expr.genes(g)
    a[[ s ]] = tf.vector(g, g.random)
  }
  a
}

# Tests motif enrichment several ways.
# Returns: a data frame of how significantly
#   different motifs were enriched, with various
#   definitions of expression enrichment.
motif.enrichment.compare = function(a, cutoff) {
  r = rbind(
    data.frame(
      method = "enriched vs depleted", cutoff = cutoff,
      enrich.test(enriched.vs.depleted(a, cutoff))),
    data.frame(
      method = "depleted vs enriched", cutoff = cutoff,
      enrich.test(enriched.vs.depleted(-a, cutoff))),

    data.frame(
      method = "enriched vs no change", cutoff = cutoff,
      enrich.test(enriched.vs.no.change(a, cutoff))),
    data.frame(
      method = "depleted vs no change", cutoff = cutoff,
      enrich.test(enriched.vs.no.change(-a, cutoff))),

    data.frame(
      method = "enriched vs random", cutoff = cutoff,
      enrich.test(enriched.versus.random(a, cutoff))),
    data.frame(
      method = "depleted vs random", cutoff = cutoff,
      enrich.test(enriched.versus.random(-a, cutoff))) )

  r[,1] = as.character(r[,1])
  r[,3] = as.character(r[,3])
  r[,4] = as.character(r[,4])
  r
}


a = NULL
for(cutoff in c(2, 1, 0.8))
  a = rbind(a,
    motif.enrichment.compare(r.sort.only.averaged, cutoff))

write.tsv(a,
  file="git/sort_paper/tf/motif/hyperg/facsEnrich/optimize1.tsv")


