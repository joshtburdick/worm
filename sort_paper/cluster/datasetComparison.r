# Clusters using various datasets.

# library("amap")
# library("ctc")

working.dir = "git/sort_paper/cluster/datasetComparison/"

source("git/utils.r")
source("git/data/name_convert.r")
source("git/unmix/seq/cluster/writeClustersTreeView.r")
# source("git/sort_paper/tf/motif/motifEnrichment2.r")
source("git/sort_paper/tf/motif/hyperg/hypergTestGroups.r")

system(paste0("mkdir -p ", working.dir, "/clustering/"))

# Reads per million (not including the HS or RNAi data.)
rpm = as.matrix(read.tsv("git/cluster/readsPerMillion.tsv"))
rpm = rpm[ , grep("HS|RNAi", colnames(rpm), invert=TRUE, value=TRUE) ]

# Log-transformed mean and max expression - used as a measure of
# whether we saw something expressed in the FACS data.
mean.expr = apply(log2(rpm+3), 1, mean) - 4

# The read ratios.
r = as.matrix(read.tsv("git/cluster/readRatios.tsv"))

# Tests for a row not being constant.
non.const.row = function(x) {
  v = var(x, na.rm=TRUE)
  !is.na(v) && (v > 0)
}

# The FACS data, including all the genes.
r.facs.orig = r[,c(1:23)]
r.facs.orig = r[ apply(r.facs.orig, 1, non.const.row) , ]

# filter 
r = r[ apply(!is.na(r), 1, all) , ]

# read Spencer data
load("data/expression/spencer.expr.Rdata")
se = as.matrix(cbind(
  spencer.expr$expr.x.1[ , c(4:33) ],
  spencer.expr$expr.y.1[ , c(4:13) ]))
rownames(se) = spencer.expr$expr.x.1$name
colnames(se)[4] = "LE.reference"   # XXX hack
se = se[ , sort(colnames(se)) ]
se = rename.gene.names(se)

# Subtracts reference samples from Spencer data
# (and removes them.)
spencer.compare.to.reference = function(se) {
  r = NULL
  for(stage in c("EE", "LE", "L2", "L3.L4", "YA")) {
    r1 = se[ , grep(stage, colnames(se)) ] -
      se[ , paste0(stage, ".reference") ]
    r = cbind(r, r1)
  }
  r = r[ , grep("reference", colnames(r), invert=TRUE) ]

  cbind(r, L1 = se[,"L1"] - se[,"LE.reference"],
    L3 = se[,"L3"] - se[,"L3.L4.reference"],
    L4 = se[,"L4"] - se[,"L3.L4.reference"])
}
se = spencer.compare.to.reference(se)

hashimshony.full.expr = read.tsv(gzfile(
  "git/data/expr/hashimshony2014/full.expr.tsv.gz"))
hy.expr = log2(1 + hashimshony.full.expr)

# only use genes shared between these
g = intersect(rownames(r), rownames(se))
r = r[g,]
se = se[g,]
se.embryonic = se[,grep("^(EE|LE)", colnames(se))]

r.facs = r[,c(1:23)]

# XXX should probably use intersection of all of these, throughout
g1 = intersect(g, rownames(hashimshony.full.expr))



# Clusters a dataset hierarchically.
cluster.dataset = function(x, name) {
  cat(paste0(name, " "))

  # do clustering
  hr = hcluster(as.matrix(x),
    method="correlation", link="complete", nbproc=2)
  cl = cutree(hr, k=300)

  write.tsv(cl, file=paste0(
    "git/sort_paper/cluster/datasetComparison/clustering_check/",
      name, ".tsv"))
}

# Clusters several datasets.
cluster.datasets = function() {
  cluster.dataset(r.facs, "facs")
  cluster.dataset(r, "facsTs")
  cluster.dataset(se, "spencer")
  cluster.dataset(se.embryonic, "spencerEmbryonic")
  cluster.dataset(hy.expr[g1,], "hy")

  if (FALSE) {

    cluster.dataset(cbind(r.facs, se), "facs.spencer")
    cluster.dataset(cbind(r, se), "facsTs.spencer")
    cluster.dataset(cbind(r.facs, se.embryonic),
      "facs.spencerEmbryonic")
    cluster.dataset(cbind(r, se.embryonic),
      "facsTs.spencerEmbryonic")


    cluster.dataset(cbind(r.facs[g1,], hy.expr[g1,]), "facs.hy")
    cluster.dataset(cbind(seEmbryonic[g1,], hy.expr[g1,]),
      "spencerEmbryonic.hy")
    cluster.dataset(cbind(se[g1,], hy.expr[g1,]), "spencer.hy")
    cluster.dataset(cbind(r.facs[g1,], seEmbryonic[g1,], hy.expr[g1,]),
      "facs.spencerEmbryonic.hy")
    cluster.dataset(cbind(r.facs[g1,], se[g1,], hy.expr[g1,]),
      "facs.spencer.hy")
  }

}

# Shuffles the numbers within each row of a dataset.
shuffle.each.row = function(r) {
  for(i in 1:nrow(r)) {
    # sample() used like this just permutes its args
    r[i,] = sample(r[i,])
  }

  r
}

# Shuffles the data, clusters it, and does enrichment analysis.
# FIXME just write out the clustering?
shuffle.cluster.enrich = function(r, output.base, n) {
  system(paste0("mkdir -p ", working.dir, "clustering/", output.base))
  motif.set.name = "hughes"

  for(iter in c(1:n)) {
    cat(iter, "\n")

    # shuffle and cluster
    cluster.dataset(shuffle.each.row(r),
      paste0(output.base, "/", iter))

    # read in clustering
    cl = as.matrix(read.tsv(paste0(working.dir, "/clustering/",
      output.base, "/", iter, ".tsv")))[,1]

    # compute enrichment, and save it
    system(paste0("mkdir -p ", working.dir, motif.set.name, "/", output.base))
    enrich = motif.enrich.many(cl, paste0("git/tf/motif/count/upstreamMotifCount/", motif.set.name, "/"))
    enrich[,,,,,"p.corr"] = p.adjust(enrich[,,,,,"p"], method="fdr")
    save(enrich, file=paste0(working.dir,
      motif.set.name, "/", output.base, "/", iter, ".Rdata"))
  }

}



# Does enrichment analysis for one clustering.
# Args:
#   motif.set.name - the name of the set of motifs
#   clustering.dir - directory containing clusterings
#   output.dir - where to write output
# Side effects: writes enrichment analysis
enrichment.analysis.1 =
    function(motif.set.name, clustering.dir, output.dir) {
  system(paste0("mkdir -p ", output.dir))

  for(clustering.file in list.files(clustering.dir)) {

    clustering.name = sub(".tsv", "", clustering.file)
    cat(paste0("\n", clustering.name, "\n"))

    # read in clustering
    cl = as.matrix(read.tsv(
      paste0(clustering.dir, "/", clustering.file)))[,1]

    enrich = motif.enrich.many(cl, paste0("git/tf/motif/count/upstreamMotifCount/", motif.set.name, "/"))
    enrich[,,,,,"p.corr"] = p.adjust(enrich[,,,,,"p"], method="fdr")
    save(enrich, file=paste0(output.dir, clustering.name, ".Rdata"))
  }
}



# cluster.datasets()

# enrichment.analysis("hughes")
# enrichment.analysis("facs", "5kb") XXX deprecated

# shuffle.cluster.enrich(r.facs, "shuffled.facs", 10)

# shuffle.cluster.labels()
if (FALSE) {
  enrichment.analysis.1("hughes",
    "git/sort_paper/cluster/datasetComparison/clustering/shuffled.labels.facs/",
    "git/sort_paper/cluster/datasetComparison/hughes/shuffled.labels.facs")
}


# shuffle and cluster
# cluster.dataset(shuffle.each.row(r.facs.orig),
#   "shuffled.rows.facs")


