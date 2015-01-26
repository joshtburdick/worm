# Looks for motifs enriched in particular clusters.

source("git/utils.r")

r = read.tsv("git/cluster/readRatios.tsv")

#motifs.1kb = read.table(
#  gzfile("git/tf/motif/motifCount/motifs_1kbUpstream.tsv.gz"),
#  sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)

# ??? not sure this is the right measure
stderr <- function(x) sqrt(var(x)/length(x))

# read in data
if (TRUE) {

#  source("git/tf/motif/motifCount/motif.counts.r")
  load("git/tf/motif/motifCount/motif.counts.Rdata")

  chip = read.table(gzfile("git/tf/chip/TF_chip_5kbUp.tsv.gz"),
    sep="\t", header=TRUE, row.names=1, check.names=FALSE,
    stringsAsFactors=FALSE)
  chip.0.5.cons = read.table(gzfile("git/tf/chip/TF_chip_5kbUp_0.5cons.tsv.gz"),
    sep="\t", header=TRUE, row.names=1, check.names=FALSE,
    stringsAsFactors=FALSE)

  # only keep genes which are in both lists
#  g = intersect(rownames(r), rownames(known.motifs))
#  r = r[g,]
#  known.motifs = known.motifs[g,]
}

# orthology data
if (TRUE) {
# instances of orthologs
ortho = read.table("git/tf/motif.ortholog.2.tsv", header=TRUE, as.is=TRUE)

# XXX check for other gene names?
# ortho = ortho[ ortho$gene %in% rownames(known.motifs), ]

# just keep cases in which there's a motif, and only keep
# the motif from the "best" homolog (by "gene.score"; arguably
# these should be normalized in some way)
ortho = ortho[ ortho$motif != "" , ]
ortho = ortho[
  order(ortho$gene,
    ortho$species,
    -ortho$related.gene.score) , ]

# XXX some of these motifs aren't included (I'm not sure why),
# so omit those cases
ortho = ortho[ ortho$motif %in% colnames(known.motifs) , ]

# ortho = ortho[ !duplicated(ortho[ , c("gene", "species"),]) , ]
# ??? just keeping unique gene names
ortho = ortho[ !duplicated(ortho[ , c("gene"),]) , ]

# tfs = intersect(rownames(known.motifs), ortho$gene)
}



if (FALSE) {

load("git/unmix/seq/cluster/hierarchical/cluster2.Rdata")
a = cutree(cluster2$hr, k=100)  # XXX arbitrary cutoff

clustering = data.frame(gene=names(a), cluster=a)
rownames(clustering) = names(a)
write.table(clustering[order(clustering$cluster),"cluster",drop=FALSE],
  file="git/unmix/seq/cluster/hierarchical/cluster1.tsv",
  sep="\t", col.names=NA)


# only use genes in both tables
g = intersect(rownames(clustering), rownames(motifs))
clustering = clustering[g,]
motifs = motifs[g,]
}



# Does many Wilcox tests.
# Args:
#   a, b - two matrices with the same number of columns
#   (NB: this will only test whether a > b)
# Returns: data frame with columns:
#   name - which column was being tested
#   W, p - results from t-test
#   a.mean, a.sd, b.mean, b.sd - motif stats in cluster
#     being tested, and background (all other clusters)
wilcox.test.many = function(a, b) {
  r = NULL

  for(j in colnames(a)) {
#    m = t.test(a[,j], b[,j], alternative="greater")
#    r = rbind(r,
#      c(name=j, m[[1]], p=m[[3]], m[[2]]))
    m = wilcox.test(a[,j], b[,j])
    r = rbind(r,
      data.frame("name"=j, "W"=m$stat, "p"=m$p.value,
        a.mean = mean(a[,j]), a.sd = sd(a[,j]),
        b.mean = mean(b[,j]), b.sd = sd(b[,j]), stringsAsFactors=FALSE))
  }

  # XXX type conversion hack
  data.frame(name = r[,"name"],
    stat = as.numeric(r[,"W"]),
    p = as.numeric(r[,"p"]),
    a.mean = as.numeric(r[,"a.mean"]),
    a.sd = as.numeric(r[,"a.sd"]),
    b.mean = as.numeric(r[,"b.mean"]),
    b.sd = as.numeric(r[,"b.sd"]))
}

# Does many t-tests, graphing the distributions.
# Args:
#   a, b - two matrices with the same number of columns
# Returns: data frame with columns:
#   name - which column was being tested
#   t, df, p - results from two-sided t-test
#   output.dir, base.name - these define the output file
t.test.many.graphing = function(a, b, output.dir, base.name) {
  system(paste("mkdir -p", output.dir))
  r = NULL

  for(j in colnames(a)) {   # XXX testing
    pos = a[,j]
    neg = b[,j]
    m = t.test(pos, neg)

    # graphing
    # XXX going by nominal p-value
    p = m[[3]]
    if (p <= 1e-5) {
      base.name = gsub(" ", "_", gsub("/", "_", base.name))
      output.file = paste(output.dir, "/", base.name, "_", j, ".pdf", sep="")
      pdf(output.file, width=2.5, height=4)
#      par(mfrow=c(1,2))
      pos.mean = mean(pos)
      neg.mean = mean(neg)
      pos.se = stderr(pos)        # 1.96 * sd(pos)
      neg.se = stderr(neg)        # 1.96 * sd(neg)
      ylim = c(min(0, pos.mean - pos.se, neg.mean - neg.se) - 0.1,
        max(0, pos.mean + pos.se, neg.mean + neg.se) + 0.1)
      mp = barplot(c(neg.mean, pos.mean),
        col=c("#4040ff", "#ff4040"), space=0.5, ylim=ylim)
      arrows(mp[1], neg.mean-neg.se, mp[1], neg.mean+neg.se,
        code=3, angle=90, lwd=2, length=0.13)
      arrows(mp[2], pos.mean-pos.se, mp[2], pos.mean+pos.se,
        code=3, angle=90, lwd=2, length=0.13)
      dev.off()
    }

    r = rbind(r,
      c(name=j, m[[1]], p=m[[3]], m[[2]]))
  }

  # XXX type conversion hack
  data.frame(name = r[,"name"],
    t = as.numeric(r[,"t"]),
    p = as.numeric(r[,"p"]),
    df = as.integer(r[,"df"]))
}

# Computes Wilxocon test for all the clusters.
# Args:
#   cluster - which cluster each gene is in
#   motif - the motifs or other predictors.
# Returns: data frame of results.
wilcox.test.all.clusters = function(cluster, motif) {
  r = NULL

  for(cl in sort(unique(cluster))) {
cat(cl, "")
    i = cluster == cl

    # skip cases in which a cluster is too small
    if (sum(i) <= 1)
      next

    z1 = (motif[i,])
    z2 = (motif[!i,])

    a = wilcox.test.many(motif[i,], motif[!i,])
#    a = t.test.many(motifs[i,], motifs[ sample(which(!i), sum(i)) , ])
    r = rbind(r, cbind(cluster = cl, a))
  }

  r
}

# Same as above, but only tests potential regulation
# by TFs associated with a cluster.
# Args:
#   expr - the expression dataset
#   cluster - which cluster each gene is in
#   motif - the motifs or other predictors
#   reg - data frame with columns "gene" and "regulator",
#     giving presumed regulatory motifs for each gene
#   min.tf.cor - only TFs at least this correlated with
#     the average of the cluster are considered 
# Returns: data frame of results.
t.test.all.clusters.filtered =
    function(expr, cluster, motif, reg, min.tf.cor) {

  r.center = compute.cluster.centers(expr, cluster, motifs)

  r = NULL

  for(cl in sort(unique(cluster))) {
cat("cluster =", cl, "\n")

    # which genes are in this cluster
    i = cluster == cl

    # which TFs are correlated with this cluster's center
    # above the cutoff
    cor.tfs = which( r.center$tf.cluster.cor[,cl] >= min.tf.cor )
    if (length(cor.tfs) == 0)
      next
cat("num. somewhat correlated TFs =", length(cor.tfs), "\n")

    # motifs associated with those TFs
    m = reg[ cor.tfs , "regulator" ]
cat("num motifs =", length(m), "\n")
    if (length(m) == 0)
      next

cat("t.test.many inputs have size", dim(motif[i,m]),
  "and", dim(motif[!i,m]), "\n")

    a = t.test.many(motif[i,m], motif[!i,m])
#    a = t.test.many(motifs[i,], motifs[ sample(which(!i), sum(i)) , ])
# print(a)

    if(nrow(a) > 0) {
      r = rbind(r, cbind(cluster=cl, a))
    }
  }

  r
}

# Writes out enrichment of only regulatory elements associated
# with TFs correlated some amount with cluster centers.
filter.correlation.cluster.enrichment =
    function(expr, motif, cluster.dir, min.tf.cor) {

  # read in the clustering
  cl = read.table(
    "git/cluster/hierarchical/hier.ts.50.clusters/clusters.tsv",
    sep="\t", header=TRUE, row.names=1, as.is=TRUE)
  cl1 = cl$cluster
  names(cl1) = rownames(cl)

  # read in the previous results
  a = read.table(paste(cluster.dir, "/", "knownMotifEnrichment_5kb.tsv", sep=""),
    sep="\t", header=TRUE, row.names=1, as.is=TRUE)

  # compute cluster centers, and correlation of motifs with same
  r.center = compute.cluster.centers(expr, cl1, known.motifs)

  # find TFs at least somewhat correlated with genes in each cluster
  r = NULL
  for(i in sort(unique(cl1))) {

    # which TFs are correlated with this cluster's center
    # above the cutoff
    cor.tfs = rownames(r.center$tf.cluster.cor)[which( r.center$tf.cluster.cor[,i] >= min.tf.cor )]
    if (length(cor.tfs) == 0)
      next
# cat("num. somewhat correlated TFs =", length(cor.tfs), "\n")
# cat(cor.tfs, "\n")
    # motifs associated with those TFs
    m = ortho[ ortho[,"gene"] %in% cor.tfs , "motif" ]
# cat("num motifs =", length(m), "\n")
    if (length(m) == 0)
      next

# cat("t.test.many inputs have size", dim(motif[i,m]),
#   "and", dim(motif[!i,m]), "\n")

    a1 = a[ a[,"cluster"]==i & (a[,"name"] %in% m) , ]

    if(nrow(a1) > 0) {
      r = rbind(r, a1)
    }

  }

  list(cl1 = cl1, r = r)
}

# Looks for enriched motifs, and writes them in a file.
# Args:
#   path - where to find "clusters.tsv" file
#   x - matrix of predictors
#   output.name - what to name the output file
# Side effects: writes a file of results
wilcox.test.clusters = function(path, x, output.name) {
  cl = read.table(paste(path, "/clusters.tsv", sep=""),
    sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
  rownames(cl) = cl$gene

  g = intersect(rownames(cl), rownames(x))
  x = x[g,]
  cl = cl[g,]

  r = wilcox.test.all.clusters(cl$cluster, x)

  r$p.fdr = p.adjust(r$p, method="fdr")
  r = r[ r$p.fdr <= 0.05 , ]
  r = r[ order(r$p.fdr) , ]

  write.table(r,
    file=paste(path, "/", output.name, ".tsv", sep=""),
    sep="\t", row.names=TRUE, col.names=NA)

  r
}



cluster.names = c(
  paste("hier.", c(50, 100, 150, 200, 250, 300), sep=""),
  paste("hier.ts.", c(50, 100, 150, 200, 250, 300), sep=""))

# using new location
cluster.dirs = paste("git/cluster/hierarchical/", cluster.names, ".clusters", sep="")
# cluster.dirs = cluster.dirs[c(10,4)]    # for testing

# Does the above test, for all clusters.
# Args:
#   motif - a set of predictors
#   output.name - name of output file
wilcox.test.all.clusters.1 = function(motif, output.name) {

  for(cluster.dir in rev(cluster.dirs)) {
    wilcox.test.clusters(cluster.dir, motif, output.name)
  }
#t.test.clusters("git/unmix/seq/cluster/hierarchical/hier",
# motif, "motifEnrichment_5kb")
#t.test.clusters("git/unmix/seq/cluster/hierarchical/hier.ts",
#  motif, "motifEnrichment_5kb")
#t.test.clusters("git/unmix/seq/cluster/WGCNA/wnet",
#  motif, "motifEnrichment_5kb")
#t.test.clusters("git/unmix/seq/cluster/WGCNA/wnet.ts",
#  motif, "motifEnrichment_5kb")

}

# Counts how many potential regulatory elements were
# enriched near genes in each cluster.
count.enriched = function(f, bh.cutoff = 0.5) {
  r = read.table(f, sep="\t", header=TRUE,
    row.names=1, stringsAsFactors=FALSE)
  r = r[ r[,"t"] > 0 & r[,"p.bh"] <= bh.cutoff, ]

  c(num = nrow(r),
    num.predictors = length(unique(r$name)),
    num.clusters = length(unique(r$cluster)))
}

###
# Does Wilcoxon tests of known motifs, de novo motifs, and ChIP
if (FALSE) {
wilcox.test.all.clusters.1(known.motifs, "MW_knownMotifEnrichment_5kb")
wilcox.test.all.clusters.1(de.novo.motifs, "MW_deNovoMotifEnrichment_5kb")
#t.test.all.clusters.1(chip, "chipEnrichment_5kb")
#t.test.all.clusters.1(known.motifs.0.5.cons, "knownMotifEnrichment_5kb_0.5cons")
#t.test.all.clusters.1(de.novo.motifs.0.5.cons, "deNovoMotifEnrichment_5kb_0.5cons")
#t.test.all.clusters.1(chip.0.5.cons, "chipEnrichment_5kb_0.5cons")
}

if (TRUE) {
  load(file="git/tf/motif/motifCount/shuffled.motif.counts.Rdata")
  wilcox.test.all.clusters.1(shuffled.known.motifs,
    "MW_shuffledKnownMotifEnrichment_5kb")
  wilcox.test.all.clusters.1(shuffled.known.motifs.0.5.cons,
    "MW_shuffledKnownMotifEnrichment_5kb_0.5cons")
}


# z = t.test.many(motif[c(1:2),], motif[c(100:120),])


# a = compute.interaction.stats("git/unmix/seq/cluster/hierarchical/hier/motifEnrichment_5kb.tsv")
# print(a)

# enriched.in.fraction()


# testing filtered-by-correlation search

cl = read.table(
  "git/cluster/hierarchical/hier.ts.200.clusters/clusters.tsv",
  sep="\t", header=TRUE, row.names=1, as.is=TRUE)
cl1 = as.character(cl$cluster)
names(cl1) = cl$gene

if (FALSE) {
r.center = compute.cluster.centers(r, cl1, known.motifs)

z = t.test.all.clusters.filtered(r, cl1, known.motifs,
  data.frame(gene=ortho[,"gene"], regulator=ortho[,"motif"]), 0.5)
z$p.corr = p.adjust(z$p, method="BH")

}

# Counts number of significant results
count.filtered.upstream = function(num.clusters, cor.cutoff) {
  r = filter.correlation.cluster.enrichment(r, known.motifs,
    paste("git/cluster/hierarchical/hier.ts.", num.clusters, ".clusters/", sep=""), cor.cutoff)
  a = r$r
  a1 = a[ a[,"p.bh"] <= 0.05,]

  c(num.clusters = num.clusters, cor.cutoff = cor.cutoff,
    num.significant = nrow(a1),
    num.motifs = length(unique(a1[,"name"])),
    num.clusters = length(unique(a1[,"cluster"])))
}


