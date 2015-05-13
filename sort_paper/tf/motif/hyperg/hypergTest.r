# Does the hypergeometric test for some clusterings.

source("git/tf/motif/enrichment/motifHyperg.r")

# XXX hokey way of getting this
orig.motif.list = {
  load("git/sort_paper/tf/allResults/motif/hier.300.clusters.Rdata")
  dimnames(enrich)[[1]]
}

# Looks for motif enrichments for one clustering.
# Args:
#   cl - the clustering
#   name - name to use for the output file
# Side effects: saves motif enrichments to that file.
enrich.test.one.clustering = function(cl, name) {
  output.dir = "git/sort_paper/tf/motif/hyperg/allResults_check/"

if (TRUE) {
  system(paste0("mkdir -p ", output.dir, "chip/"))
  enrich = enrich.test.many.motifs(
    "git/tf/chip/count/upstreamChipCount/", cl)
  save(enrich, file=
    paste0(output.dir, "chip/", name, ".Rdata"))
}
if (TRUE) {
  system(paste0("mkdir -p ", output.dir, "5kb/"))
  enrich = enrich.test.many.motifs(
    "git/tf/motif/count/upstreamMotifCount/5kb/", cl, orig.motif.list)
  save(enrich, file=
    paste0(output.dir, "5kb/", name, ".Rdata"))

  system(paste0("mkdir -p ", output.dir, "hughes/"))
  system("mkdir -p git/sort_paper/tf/motif/hyperg/allResults/hughes/")
  enrich = enrich.test.many.motifs(
    "git/tf/motif/count/upstreamMotifCount/hughes_20141202/", cl)
  save(enrich, file=
    paste0(output.dir, "hughes/", name, ".Rdata"))
}
if (FALSE) {
  system(paste0("mkdir -p ", output.dir, "jolma2013_shuffled/"))
  enrich = enrich.test.many.motifs(
    "git/tf/motif/count/upstreamMotifCount/jolma2013_shuffled/", cl, orig.motif.list)
  save(enrich, file=
    paste0(output.dir, "jolma2013_shuffled/", name, ".Rdata"))
}
}

if (TRUE) {
  # go through the directory of clusterings
  for(clustering.name in list.files("git/cluster/hierarchical")) {
    write.status(clustering.name)

    cl = as.matrix(read.tsv(
      paste0("git/cluster/hierarchical/", clustering.name, "/clusters.tsv")))[,2]
    cl = gsub(" ", "", cl)

    enrich.test.one.clustering(cl, clustering.name)
  }

#  enrich = enrich.test.many.motifs(
#    "git/tf/motif/count/upstreamMotifCount/hughes_20141202/", cl)

#  save(enrich, file=
#    "git/sort_paper/tf/motif/hyperg/hughes_20141202/facs.Rdata")
}

if (FALSE) {

  # one way of shuffling, by just shuffling the cluster labels
  cl.shuffled = sample(cl)
  names(cl.shuffled) = names(cl)

  shuffled.enrich = enrich.test.many.motifs(
     "git/tf/motif/count/upstreamMotifCount/hughes/", cl.shuffled)
  save(shuffled.enrich, file=
    "git/sort_paper/tf/motif/hyperg/shuffled.enrich.Rdata")

  # Shuffling the expression dataset, then clustering
  # XXX this is awkward
  cl.shuffled.expr = {
    r = read.tsv(
    "git/sort_paper/cluster/datasetComparison/clustering/shuffled.rows.facs.tsv")
    r1 = as.character(r[,1])
    names(r1) = rownames(r)
    r1
  }

shuffled.expr.enrich = enrich.test.many.motifs(
  "git/tf/motif/count/upstreamMotifCount/hughes/", cl.shuffled.expr)
save(shuffled.expr.enrich, file=
  "git/sort_paper/tf/motif/hyperg/shuffled.expr.enrich.Rdata")
}

# Runs on one directory.
enrich.test.file = function(f) {
  output.dir = "git/sort_paper/tf/motif/hyperg/hughes_20141202/"
  system(paste0("mkdir -p ", output.dir))

  cl = as.matrix(read.tsv(
    paste0("git/sort_paper/cluster/datasetComparison/clustering/",
      f)))[,1]
  # XXX
  cl = { 
    cl.names = names(cl)
    cl = as.character(cl)
    names(cl) = cl.names
    cl
  }
  name = sub(".tsv", "", f)

  enrich.test.one.clustering(cl, name)

#  enrich = enrich.test.many.motifs(
#    "git/tf/motif/count/upstreamMotifCount/hughes_20141202/", cl)
#
#  save(enrich, file=paste0(output.dir, "/", name, ".Rdata"))
}

# Runs on one directory of the original clusterings.
enrich.test.clustering = function() {
  output.dir = "git/sort_paper/tf/motif/hyperg/hughes_20141202_facs_clustering/"
  system(paste0("mkdir -p ", output.dir))

  for(clustering.name in list.files("git/cluster/hierarchical")) {
    cat(paste0(clustering.name, "\n"))
    cl = as.matrix(read.tsv(
      paste0("git/cluster/hierarchical/", clustering.name, "/clusters.tsv")))[,1]
    cl = gsub(" ", "", cl)

    enrich = enrich.test.many.motifs(
      "git/tf/motif/count/upstreamMotifCount/hughes_20141202/", cl)

    save(enrich, file=paste0(output.dir, "/", clustering.name, ".Rdata"))
  }
}

if (FALSE) {
  enrich.test.file("hy.tsv")
  enrich.test.file("spencer.tsv")
  enrich.test.file("spencerEmbryonic.tsv") 
  enrich.test.file("shuffled.rows.facs.tsv")
#  enrich.test.clustering()
}


