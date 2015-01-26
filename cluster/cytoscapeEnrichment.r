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
if (FALSE) {
ortho = ortho[
  order(ortho$gene,
    ortho$species,
    -ortho$related.gene.score) , ]
ortho = ortho[ !duplicated(ortho[ , c("motif"),]) , ]
}

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

# chip data and de novo motifs (omitted for now)
if (FALSE) {
  r1 = read.tsv(paste(path, "chipEnrichment_5kb.tsv", sep="/"))
  r1$type = "ChIP"
  r1$gene = chip.gene.name[ r1$name , "gene" ]

  # "fold enrichment", in this case, is just difference in means
  r1$log2.fold.enrich = r1$a.mean - r1$b.mean

  r1 = r1[ , c("name", "type", "cluster", "p.bh", "gene", "a.mean","b.mean","log2.fold.enrich")]

#  r2 = read.tsv(paste(path, "deNovoMotifEnrichment_5kb.tsv", sep="/")) 
#  r2$type = "deNovoMotif"
#  r2$gene = ""
#  r2 = r2[ , c("name", "type", "cluster", "p.bh", "gene")]
}

  r3 = read.tsv(paste(path, "uniqueKnownMotifEnrichment_5kb_0.5cons.tsv", sep="/"))
  r3$type = "knownMotif"
  r3$log2.enrich = log2(r3$enrichment)

  # add one row for each gene with a given motif
  # (since there are more motifs than genes)
  # XXX this is slow, but, um, hopefully correct
if (FALSE) {
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
  r3a$log2.fold.enrich = log2( r3a$enrichment )

cat("\n")
  r3a = r3a[ , c("motif", "type", "cluster", "p.fdr", "gene","log2.fold.enrich") ]
  colnames(r3a)[1] = "name"
print(dim(r3))
print(dim(r3a))
}

  r = rbind(r3)       # XXX temporary
  r$edge_name = r$name

  r = r[ r$p.fdr <= 0.05 , ]

  # add on correlation of gene with each TF
  clusters = read.tsv(paste(path, "/clusters.tsv", sep=""))
  cc = compute.cluster.centers(expr.ratio, clusters$cl, unique(r$gene))
  i = r$gene %in% rownames(cc$tf.cluster.cor)
  r$tf.correlation = NA
  r$tf.correlation[i] = cc$tf.cluster.cor[ cbind(r$gene[i], r$cluster[i]) ]

  # ??? filter by correlation?
#  r = r[ !is.na(r$tf.correlation) , ]
#  r = r[ abs(r$tf.correlation) >= 0.5 , ]

  # filter by "enrichment"
#  r = r[ abs(r$log2.fold.enrich) >= 0.5 , ]
  r = r[ r$log2.fold.enrich >= 0.5 , ]

  write.tsv(r, paste(path, "cytoscape.tsv", sep="/"))
}

write.motif.enrichment.cytoscape.file = function(path) {
  r3 = read.table(
    paste(path, "uniqueKnownMotifEnrichment_5kb_0.5cons.tsv", sep="/"),
    sep="\t", header=TRUE, as.is=TRUE)
  r3$type = "knownMotif"
  r3$log2.fold.enrich = log2(r3$enrichment)
  r3$log2.fold.enrich[ r3$log2.fold.enrich <= -100 ] = -100
  
  r3 = r3[ abs(r3$log2.fold.enrich) >= log2(2) , ]
  r3 = r3[ , c("motif", "cluster", "p.fdr", "type", "log2.fold.enrich") ]
  colnames(r3) = c("a", "b", "note", "type", "strength")
  # scale this so that an eight-fold enrichment is + 1
  r3$strength = r3$strength / 3

#  r.ortho = ortho[ ortho$motif %in% unique(r3$motif) , ]
  # get orthologs
  r.ortho = ortho
  r.ortho = data.frame(a=r.ortho$gene, b=r.ortho$motif,
    note="", type="orthology", strength=0)
  r.ortho = r.ortho[ r.ortho$b %in% r3$a , ]
  r.ortho = unique(r.ortho)
print(dim(r.ortho))

  # add on correlation of gene with each TF
  clusters = read.tsv(paste(path, "/clusters.tsv", sep=""))
  cc = compute.cluster.centers(expr.ratio, clusters$cl, unique(r.ortho$a))
#  i = r.ortho$a %in% rownames(cc$tf.cluster.cor)
print(cc$tf.cluster.cor[1:5,1:5])
  tf.cor = cc$tf.cluster.cor
  r.cor = NULL
  for(j in colnames(tf.cor)) {
    r.cor = rbind(r.cor,
      data.frame(a=rownames(tf.cor), b=as.character(j), note="",
        type="tfCorrelation",
        strength = tf.cor[,as.character(j)]))
  }

#  r.cor = data.frame(a=r3$a[i], b=r3$b[i], note="",
#    type="tfCorrelation",
#    strength=(1/0.8) * cc$tf.cluster.cor[ cbind(r.ortho$a[i], as.character(r.ortho$b[i])) ])
  r.cor = r.cor[ abs(r.cor$strength) >= 0.6 , ]

  r = rbind(r3, r.ortho, r.cor)

  write.tsv(r, paste(path, "cytoscape.tsv", sep="/"))
}

# Writes enrichment just for motifs.
if (TRUE) {
  write.motif.enrichment.cytoscape.file("git/cluster/hierarchical/hier.250.clusters/")
}

if (FALSE) {
write.enrichment.cytoscape.file("git/cluster/hierarchical/hier.200.clusters/")
}


