# Constructs a network based on things (motifs and ChIP signal)
# upstream of genes in clusters, filtered by correlation of a TF
# with clusters.

source("git/utils.r")

# the expression data, as ratios
expr.ratio = read.tsv("git/cluster/readRatios.tsv")

# the motif counts
load("git/tf/motif/motifCount/motif.counts.Rdata")

# orthology info
ortho = read.table("git/tf/motif.ortholog.2.tsv", header=TRUE, as.is=TRUE)
# omitting this for now (I may just not be parsing it right)
ortho = ortho[ ortho$method != "InParanoid" , ]
# XXX some of these motifs aren't included (I'm not sure why),
# so omit those cases
ortho = ortho[ ortho$motif %in% colnames(known.motifs) , ]
# just keep cases in which there's a motif, and only keep
# the "best" ortholog for each gene (by "gene.score";
# arguably these should be normalized in some way)
ortho = ortho[
  order(ortho$gene,
    ortho$species,
    -ortho$related.gene.score) , ]
ortho = ortho[ !duplicated(ortho[ , c("gene"),]) , ]

# the ChIP data
chip = read.table(gzfile("git/tf/chip/TF_chip_5kbUp.tsv.gz"),
  sep="\t", header=TRUE, row.names=1, check.names=FALSE,
  stringsAsFactors=FALSE)
chip.0.5.cons = read.table(gzfile("git/tf/chip/TF_chip_5kbUp_0.5cons.tsv.gz"),
  sep="\t", header=TRUE, row.names=1, check.names=FALSE,
  stringsAsFactors=FALSE)
chip.gene.name = read.tsv("git/tf/chip/chip.gene.name.tsv")


# Computes center of each cluster.
# Args:
#   r - the expression data
#   cl - the clustering (as a vector of integers, with genes as names)
#   tf - names of the genes to consider TFs
# Returns: list with
#   r.center - center of each cluster (for now, using the mean)
#   tf.cluster.cor - correlation of each TF wich each cluster center
compute.cluster.centers = function(r, cl, tf) {
  cluster.names = sort(unique(cl))
  num.clusters = max(cl)

  # first, compute cluster centers
  r.center = matrix(nrow=num.clusters, ncol=ncol(r))
  rownames(r.center) = cluster.names
  colnames(r.center) = colnames(r)

  for (i in cluster.names)
    r.center[i,] = apply(r[cl==i,], 2,
      function(x) mean(x, na.rm=TRUE))

  # then, correlation of each TF with each cluster
  g = intersect(rownames(r), tf)
  tf.expr = r[g,]
  tf.cluster.cor = cor(t(r[g,]), t(r.center))

  list(r.center = r.center, tf.cluster.cor = tf.cluster.cor)
}

# Reads in the various enrichments, and writes out a file
# which Cytoscape can read.
# Args:
#   path - where the various enrichment files are
# Side effects: writes a Cytoscape-readable file in that directory.
write.enrichment.cytoscape.file = function(path) {

if (FALSE) {
  r1 = read.tsv(paste(path, "chipEnrichment_5kb.tsv", sep="/"))
  r1$type = "ChIP"
  r1$gene = chip.gene.name[ r1$name , "gene" ]

  # "fold enrichment", in this case, is just difference in means
  r1$log2.fold.enrich = r1$a.mean - r1$b.mean

  r1 = r1[ , c("name", "type", "cluster", "p.bh", "gene", "a.mean","b.mean","log2.fold.enrich")]
}
#  r2 = read.tsv(paste(path, "deNovoMotifEnrichment_5kb.tsv", sep="/")) 
#  r2$type = "deNovoMotif"
#  r2$gene = ""
#  r2 = r2[ , c("name", "type", "cluster", "p.bh", "gene")]

  r3 = read.tsv(paste(path, "chisq_knownMotifEnrichment_5kb.tsv", sep="/"))
  r3$type = "knownMotif"

  # add one row for each gene with a given motif
  # (since there are more motifs than genes)
  # XXX this is slow, but, um, hopefully correct
  r3a = NULL
cat("\n")
  for (i in 1:nrow(ortho)) {
cat(backspace.string, i)
    a = r3[ r3$motif == ortho[i,"motif"] , ]
    if (nrow(a) > 0) {
      a$gene = ortho[ i, "gene" ]
      r3a = rbind(r3a, a)
    }
  }

  # convert means to "fold enrichment"
#  r3a$log2.fold.enrich = log2(r3a$a.mean) - log2(r3a$b.mean)
#  r3a$log2.fold.enrich[ r3a$a.mean == 0 | r3a$b.mean == 0 ] = 0
  r3a$log2.fold.enrich = log2( r3a$enrichment )

cat("\n")
  r3a = r3a[ , c("motif", "type", "cluster", "p.fdr", "gene","log2.fold.enrich") ]
  colnames(r3a)[1] = "name"
print(dim(r3))
print(dim(r3a))

  r = rbind(r3a)       # XXX temporary
  r$edge_name = r$name

  # XXX this doesn't seem to help much
  r = r[ r$p.fdr <= 0.05 , ]

  # add on correlation of gene with each TF
  clusters = read.tsv(paste(path, "/clusters.tsv", sep=""))
  cc = compute.cluster.centers(expr.ratio, clusters$cl, unique(r$gene))
  i = r$gene %in% rownames(cc$tf.cluster.cor)
  r$tf.correlation = NA
  r$tf.correlation[i] = cc$tf.cluster.cor[ cbind(r$gene[i], r$cluster[i]) ]

  # ??? filter by correlation?
  r = r[ !is.na(r$tf.correlation) , ]
#  r = r[ abs(r$tf.correlation) >= 0.5 , ]

  # filter by "enrichment"
#  r = r[ abs(r$log2.fold.enrich) >= 0.5 , ]
  r = r[ r$log2.fold.enrich >= 0.5 , ]

  write.tsv(r, paste(path, "cytoscape.tsv", sep="/"))
}

write.motif.enrichment.cytoscape.file = function(path) {
  r3 = read.table(
    paste(path, "uniqueKnownMotifEnrichment_5kb.tsv", sep="/"),
    sep="\t", header=TRUE, as.is=TRUE)
  r3$type = "knownMotif"
  r3$log2.fold.enrich = log2(r3$enrichment)
  r3$log2.fold.enrich[ r3$log2.fold.enrich <= -100 ] = -100
  r3 = r3[ r3$log2.fold.enrich >= log2(2) , ]

  r = r3

  write.tsv(r, paste(path, "cytoscape.tsv", sep="/"))
}

# Writes enrichment just for motifs.
if (TRUE) {
  write.motif.enrichment.cytoscape.file("git/cluster/hierarchical/hier.100.clusters/")
}

if (FALSE) {
write.enrichment.cytoscape.file("git/cluster/hierarchical/hier.ts.200.clusters/")
}


