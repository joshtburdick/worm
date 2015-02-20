# Quick look at predictions for sets of genes with
# particular binding sites.

source("git/sort_paper/unmix/pseudoinverseEnrichment.r")

# XXX hokey way of getting this
orig.motif.list = {
  load("git/sort_paper/tf/allResults/motif/hier.300.clusters.Rdata")
  dimnames(enrich)[[1]]
}

# Gets particular motifs.
# Args:
#   motif.dir - directory containing the motifs
#   motifs - the motifs to get
get.motif.counts = function(motif.dir, motifs) {

  get.counts = function(m) {
    motif.count = NULL
    f = paste0(motif.dir, "/", m, ".Rdata")
    if (file.exists(f))
      load(f)
    motif.count[,"1","0.9","35"]
  }

  m1 = get.counts(motifs[1])
  m = matrix(nrow = length(m1), ncol=length(motifs))
  dimnames(m) = list(gene=names(m1), motif=motifs)

  m[,1] = m1
  for(j in 2:length(motifs))
    m[,j] = get.counts(motifs[j])

  m
}

motif.counts = get.motif.counts(
  "git/tf/motif/count/upstreamMotifCount/5kb/",
  orig.motif.list)

g = intersect(rownames(x.pseudoinverse), rownames(motif.counts))
motif.counts = motif.counts[g,]
xp1 = x.pseudoinverse[g,]

x.with.motif = t(motif.counts > 0) %*% xp1

image(t(x.with.motif))




