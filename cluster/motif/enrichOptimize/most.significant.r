# Gets the cutoffs which were most significant.

source("git/utils.r")

motif.info.dir = "git/cluster/motif/distAndConservation/5kb/"
clustering.name = "hier.300.clusters"

# table for renaming motifs
motif.names = read.table("data/tf/meme/motifList.tsv",
  sep="\t", header=TRUE, as.is=TRUE)
motif.names[ motif.names$name == "" , "name" ] = motif.names[ motif.names$name=="", "id"]
motif.id.to.name = motif.names$name
names(motif.id.to.name) = motif.names$id

# read in clustering
cl1 = read.tsv(paste0("git/cluster/hierarchical/", clustering.name, "/clusters.tsv"))
colnames(cl1) = c("gene", "set")
cl = cl1$set
names(cl) = cl1$gene

# For each motif and cluster, finds the most significant
# p-value, and the settings with which that p-value
# was found.
# XXX deprecated -- this code is fairly opaque.
motif.enrich.most.sig = function(clustering.name) {

  # read in clustering
  cl1 = read.tsv(paste0("git/cluster/hierarchical/", clustering.name, "/clusters.tsv"))
  colnames(cl1) = c("gene", "set")
  cl = cl1$set
  names(cl) = cl1$gene

  # enriched motifs per-cluster
  load(paste0("git/cluster/motif/enrichOptimize/cutoff.optimize/", clustering.name, ".Rdata"))

  # get ordering of which things were significant
  cm.min = apply(p.corr, c(1,2), min)
#  cm.sig.ranking = apply(p.corr, c(1,2), order)

  m = apply(p.corr, c(1,2), which.min)
  list(cm.min = apply(p.corr, c(1,2), min),
    upstream.dist.kb = m %% 3 + 1,
    conservation = (m %/% 3) %% 4 + 1,
    motif.score = (m %/% 12) %% 3 + 1)
}

# Potentially sloooooooooooow function which finds, for each gene,
# the k most highly enriched motifs upstream of it.
# Args:
#   cl - a clustering (as a numeric vector, named by gene)
#   most.sig.stats - statistics about which motifs were
#     significantly enriched (see motif.enrich.most.sig())
# Returns: data frame with columns
#   gene - which gene
#   motif - which motif
#   enrich.p.corr - the FDR-corrected enrichment p-value
# Note that the motif in question may be present upstream of the
# gene in question, outside of the parameters with which the
# enrichment was found "significant".
per.gene.motifs = function(cl, most.sig.stats) {
  gene.motif = NULL

  for(motif in rownames(most.sig.stats$cm.min)) {
    write.status(motif)

    r = read.table(paste0(motif.info.dir, "/", motif, "_upstreamMotifCons.tsv.gz"),
      sep="\t", as.is=TRUE)
    # XXX should have included column names, and maybe upstream distance also
    colnames(r) =
      c("region.chr", "region.a", "region.b",
      "gene", "score", "strand",
      "motif", "motif.a", "motif.b", "motif.id", "motif.score",
      "motif.strand", "motif.cons")
    r$upstream.dist = ifelse(r$strand=="+",
      r$motif.a - r$region.b, r$region.a - r$motif.a)

    # tack on the clustering, and the significance 
    r$cl = as.character( cl[r$gene] )
    r = r[ !is.na(r$cl) , ]
    r$enrich.p.corr = most.sig.stats$cm.min[ motif, r$cl ]

    # filter these such that we only include the motifs in the
    # groups which were enriched
    r = r[ r$upstream.dist >= most.sig.stats$upstream.dist.kb[motif,r$cl] * -1000 , ]
    r = r[ r$motif.cons >= c(0, 0.5, 0.7, 0.9)[most.sig.stats$conservation[motif,r$cl] ] , ]
    r = r[ r$motif.score >= c(30, 35, 40)[most.sig.stats$motif.score[motif,r$cl] ] , ]

    # filter by significance, and uniquify
    r = unique( r[ r$enrich.p.corr <= 0.05 , c("gene","cl","enrich.p.corr") ] )
    r = cbind(motif=motif, r, stringsAsFactors=FALSE)
    gene.motif = rbind(gene.motif, r)
  }

  gene.motif
}

# For each gene pair, finds the k most significant motifs.
find.most.significant = function(r, k) {
  f = function(x) {
    i = order(x$enrich.p.corr)
    if (length(i) >= k)
      x = x[ i[1:k] , ]
    else
      x = x[ i , ]
  }

  by(r, r$gene, f)
}

# Finds the most significant "settings" for motif enrichment, for each motif.
most.significant.cutoffs = function(most.sig.stats) { 
  j = apply(most.sig.stats$cm.min, 1, which.min)
  n = length(j)
  data.frame(motif = rownames(most.sig.stats$cm.min),
    upstream.dist.kb = most.sig.stats$upstream.dist.kb[ cbind(1:n,j) ],
    conservation = c(0,0.5,0.7,0.9)[ most.sig.stats$conservation[ cbind(1:n,j) ] ],
    motif.score = c(30,35,40)[ most.sig.stats$motif.score[ cbind(1:n,j) ] ] )
}
write.tsv(most.significant.cutoffs(most.sig.stats),
  "git/cluster/motif/enrichOptimize/mostSignificantCutoffs.hier.300.tsv")

# Summarizes the motifs.
significant.motif.summary = function(r) {
  f = function(a) {
    lp = ifelse(a$enrich.p.corr <= 1e-100, 100, -log10(a$enrich.p.corr))
    paste(paste0(motif.id.to.name[a$motif], "(", trunc(lp), ")"), collapse = " ")
  }

  sapply(r, f)
}

clustering.name = "hier.300.clusters"

# read in clustering
cl1 = read.tsv(paste0("git/cluster/hierarchical/", clustering.name, "/clusters.tsv"))
colnames(cl1) = c("gene", "set")
cl = cl1$set
names(cl) = cl1$gene

most.sig.stats = motif.enrich.most.sig(clustering.name)
write.tsv(most.significant.cutoffs(most.sig.stats),
  "git/cluster/motif/enrichOptimize/mostSignificantCutoffs.hier.300.tsv")

motifs.and.genes = per.gene.motifs(cl, most.sig.stats)

motifs.per.gene = significant.motif.summary(find.most.significant(motifs.and.genes, 5))

# write.tsv(data.frame(motif=motifs.per.gene), "~/foo.tsv")


