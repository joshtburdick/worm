# Looks for motifs enriched in particular clusters.

library(rpart)   # deprecated

#motifs = read.table(
#  gzfile("git/tf/motif/motifCount/motifs_1kbUpstream.tsv.gz"),
#  sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)

# ??? not sure this is the right measure
stderr <- function(x) sqrt(var(x)/length(x))

# read in data
if (TRUE) {

#  source("git/tf/motif/motifCount/motif.counts.r")
  load("git/tf/motif/motifCount/motif.counts.Rdata")

  chip = read.table(gzfile("git/tf/chip/TF_chip_5kbUp.tsv.gz"),
    sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
}








# Looks for enriched motifs, and writes them in a file.
# Args:
#   path - where to find "clusters.tsv" file
#   x - matrix of predictors
#   output.name - what to name the output file
# Side effects: writes a file of results
t.test.clusters = function(path, x, output.name) {
  cl = read.table(paste(path, "/clusters.tsv", sep=""),
    sep="\t", header=TRUE, row.names=1, stringsAsFactors=TRUE)
  rownames(cl) = cl$gene

  g = intersect(rownames(cl), rownames(x))
  x = x[g,]
  cl = cl[g,]

  r = t.test.all.clusters(cl$cluster, x)
  r$p.bh = p.adjust(r$p, method="hochberg")
  r = r[ r$p.bh <= 0.5 & r$t > 0 , ]
#  r = r[ r$p <= 0.01 & r$t > 0 , ]
  r = r[ order(r$p.bh) , ]

  write.table(r,
    file=paste(path, "/", output.name, ".tsv", sep=""),
    sep="\t", row.names=TRUE, col.names=NA)

  r
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

# one sort of regression...
# cat("doing regression\n")
# a = cbind(cluster = factor(clustering$cluster), motifs)
# f = paste("cluster ~ ", paste(colnames(motifs), collapse="+"))
# r = rpart(formula(f), data=a)

# Does many t-tests.
# Args:
#   a, b - two matrices with the same number of columns
# Returns: data frame with columns:
#   name - which column was being tested
#   t, df, p - results from two-sided t-test
t.test.many = function(a, b) {
  r = NULL

  for(j in colnames(a)) {
    m = t.test(a[,j], b[,j])
    r = rbind(r,
      c(name=j, m[[1]], p=m[[3]], m[[2]]))
  }

  # XXX type conversion hack
  data.frame(name = r[,"name"],
    t = as.numeric(r[,"t"]),
    p = as.numeric(r[,"p"]),
    df = as.integer(r[,"df"]))
}

# Same as above, but also graphs the distributions.
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

# Tests for motifs enriched in particular sorted fractions.
enriched.in.fraction.1 = function(log.enrich, motif, output.dir) {
  cutoff = log2(3)

  g1 = intersect(rownames(log.enrich), rownames(motif))
  log.enrich = log.enrich[g1,]
  motif = motif[g1,]
  a = NULL
  for(s in colnames(log.enrich)) {
    cat(s, "")
    r = t.test.many(motif[ log.enrich[,s] >= cutoff , ],
      motif[ log.enrich[,s] <= -cutoff , ])
#    r = t.test.many.graphing(motif[ log.enrich[,s] >= cutoff , ],
#      motif[ log.enrich[,s] <= -cutoff , ], output.dir, s)
    r$p.bh = p.adjust(r$p, method="hochberg")
    r = r[ r$p.bh <= 0.5 & r$t > 0 , ]
    r = r[ order(r$p.bh) , ]
    if (nrow(r) > 0)
      r = data.frame(sample=s, r)
    a = rbind(a, r)
  }

  a
}

# Does that test for all sorted fractions.
enriched.in.fraction = function() {
  output.path = "git/unmix/seq/cluster/fraction_enrich"
  system(paste("mkdir -p", output.path))

  r = as.matrix(read.table("git/unmix/seq/cluster/readsNormalized.tsv",
    header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

  write.table(enriched.in.fraction.1(r[,c(1:25)], motif,
    paste(output.path, "/motif_5kb", sep="")),
    file=paste(output.path, "/motif_5kb.tsv", sep=""),
    sep="\t", row.names=TRUE, col.names=NA)
  write.table(enriched.in.fraction.1(r[,c(1:25)], chip,
    paste(output.path, "/chip_5kb", sep="")),
    file=paste(output.path, "/chip_5kb.tsv", sep=""),
    sep="\t", row.names=TRUE, col.names=NA)
}

# Computes t-tests for all the clusters.
# Args:
#   cluster - which cluster each gene is in
#   motif - the motifs or other predictors.
# Returns: data frame of results.
t.test.all.clusters = function(cluster, motif) {
  r = NULL

  for(cl in sort(unique(cluster))) {
cat(cl, "")
    i = cluster == cl

  z1 = (motif[i,])
  z2 = (motif[!i,])

    a = t.test.many(motif[i,], motif[!i,])
#    a = t.test.many(motifs[i,], motifs[ sample(which(!i), sum(i)) , ])
    r = rbind(r, cbind(cluster = cl, a))
  }

  r
}

cluster.dirs = paste("git/unmix/seq/cluster/",
  c("hierarchical/hier.50clusters/",
    "hierarchical/hier.100clusters/",
    "hierarchical/hier.200clusters/",
    "hierarchical/hier.ts.50clusters/",
    "hierarchical/hier.ts.100clusters/",
    "hierarchical/hier.ts.200clusters/",
    "WGCNA/wnet_pow11_mch0.1/",
    "WGCNA/wnet_pow11_mch0.2/",
    "WGCNA/wnet.ts_pow11_mch0.1/",
    "WGCNA/wnet.ts_pow11_mch0.2/"), sep="")

# using new location
cluster.dirs = c(paste("git/cluster/hierarchical/hier.",
    c(50, 100, 150, 200, 250, 300), ".clusters", sep=""),
  paste("git/cluster/hierarchical/hier.ts.",
    c(50, 100, 150, 200, 250, 300), ".clusters", sep=""))


# Does the above test, for all clusters.
# Args:
#   motif - a set of predictors
#   output.name - name of output file
t.test.all.clusters.1 = function(motif, output.name) {

  for(cluster.dir in cluster.dirs) {
    t.test.clusters(cluster.dir, motif, output.name)
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

# Computes statistics about some clustering.
compute.interaction.stats = function(f, bh.cutoff = 0.5) {
  r = read.table(f, sep="\t", header=TRUE,
    row.names=1, stringsAsFactors=FALSE)
  r = r[ r[,"p.bh"] <= bh.cutoff , ]

  c(num = nrow(r),
    num.predictors = length(unique(r$name)),
    num.clusters = length(unique(r$cluster)))
}


if (FALSE) {
t.test.all.clusters.1(known.motifs, "knownMotifEnrichment_5kb")
t.test.all.clusters.1(de.novo.motifs, "deNovoMotifEnrichment_5kb")
t.test.all.clusters.1(chip, "chipEnrichment_5kb")

# Computes stats about enrichment for all the clusterings.
cluster.stats = function() {
  r = NULL
  for(cl in cluster.dirs[1:4]) {
    r = rbind(r,
      c(cl,
        compute.interaction.stats(
          paste(cl, "knownMotifEnrichment_5kb.tsv", sep="/"), 0.05),
        compute.interaction.stats(
          paste(cl, "deNovoMotifEnrichment_5kb.tsv", sep="/"), 0.05)))
  }
  colnames(r) = c("clustering",
    "known.motifs.pairs", "known.motifs", "known.motifs.clusters",
    "de.novo.motifs.pairs", "de.novo.motifs", "de.novo.motifs.clusters")
  r
}
write.table(cluster.stats(), file="git/cluster/clusterStats.tsv",
  sep="\t", row.names=TRUE, col.names=NA)
}

# z = t.test.many(motif[c(1:2),], motif[c(100:120),])


# a = compute.interaction.stats("git/unmix/seq/cluster/hierarchical/hier/motifEnrichment_5kb.tsv")
# print(a)

# enriched.in.fraction()

