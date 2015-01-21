# Does the motif-enrichment hypergeometric test described in
# Barash et al (2001).

source("git/utils.r")

# Looks for significantly-enriched motifs, using the
# hypergeometric test mentioned above.
# Args (which will be recycled according to R's rules):
#   m.cluster - the number of genes with the motif in each cluster
#   g.cluster - total number of genes in each cluster
#   m.total - the total number of genes with the motif
#   g.total - the total number of genes
# Returns: a vector of (uncorrected) hypergeometric p-values.
motif.enrich.hyperg =
    function(m.cluster, g.cluster, m.total, g.total) {
  phyper(m.cluster - 1, m.total, g.total - m.total, g.cluster,
    lower.tail=FALSE)
}

# The same thing, but implemented naively (as in section 2.3 of
# that paper; used for testing.)
motif.enrich.hyperg.naive =
    function(m.cluster, g.cluster, m.total, g.total) {
  p = NA
  s = 0
  if (m.cluster[i] <= m.total[i]) {   # XXX shouldn't happen
    for(k1 in m.cluster[i]:m.total[i]) {
      s = s + dhyper(k1, m.total[i],
        g.total - m.total[i], g.cluster[i])
    }
    p = s
  }
  p
}
# sample call: r = motif.enrich.hyperg(10, 100, 30, 1000)

# Tests for enrichment of many motifs.
# Args:
#   motif.dir - directory containing motif counts
#   cl - a clustering of genes
#   motifs - optional subset of motifs to include
# Returns: a data frame of enrichments
enrich.test.many.motifs = function(motif.dir, cl, motifs=NULL) {
  cluster.names = sort(unique(as.character(cl)))
  motif.files = list.files(motif.dir)    # [1:20]
  motif.names = sub(".Rdata$", "", motif.files)

  # possibly filter list of motifs to include
  if (!is.null(motifs)) {
    motif.names = intersect(motif.names, motifs)
  }
  motif.files = intersect(motif.files, paste0(motif.names, ".Rdata"))

  # get a file of counts, and various dimensions
  load(paste0(motif.dir, "/", motif.files[1]))

  g = intersect(names(cl), dimnames(motif.count)[[1]])
  cl = cl[g]
  g.cluster = c(by(cl, cl, length))
  g.cluster = g.cluster[ cluster.names ]
  g.total = length(g)

  # object to hold results
  a = list(motif = motif.names,
    group = cluster.names,
    stat=c("m.cluster", "g.cluster", "m.total", "p", "p.corr"))
  r = array.from.dimnames(c(a, dimnames(motif.count)[2:4]))

  # loop through the motif files
  if (TRUE) {
    for(f in motif.files) {
      write.status(f)
      load(paste0(motif.dir, "/", f))
      motif.name = sub(".Rdata$", "", f)

      # restrict to the genes clustered
      # XXX also using an arbitrary cutoff for "motif is present"
      motif.count = motif.count[g,,,] > 0

      sum.by.cluster = function(x) {
        c(by(x, cl, sum))
      }
      m.cluster = apply(motif.count, c(2:4), sum.by.cluster)
      m.total = apply(motif.count, c(2:4), sum)
      r1 = motif.enrich.hyperg(m.cluster, g.cluster,
        rep(m.total, each=length(g.cluster)), g.total)

      r[motif.name,,"m.cluster",,,] = m.cluster
      r[motif.name,,"g.cluster",,,] = g.cluster
      r[motif.name,,"m.total",,,] = m.total
      r[motif.name,,"p",,,] = r1
    }
  }

  r[,,"p.corr",,,] = p.adjust(r[,,"p",,,], method="fdr")
  r
}

