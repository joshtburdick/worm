# Does the hypergeometric test for some clusterings.

source("git/tf/motif/enrichment/motifHyperg.r")

# XXX hokey way of getting this
orig.motif.list = {
  load("git/sort_paper/tf/allResults/motif/hier.300.clusters.Rdata")
  dimnames(enrich)[[1]]
}

# Does enrichment tests for one clustering, and one set of motifs.
hyperg.test.1 = function(cl, motif.set.name, clustering.name) {
  motif.count.base =
    "/home/jburdick/gcb/git/sort_paper/tf/motif/upstreamCount/"
  output.dir = "git/sort_paper/tf/motif/hyperg/allResults2/"
  motif.count.dir = paste0(motif.count.base, "/", motif.set.name, "/")

  system(paste0("mkdir -p ", output.dir, "/", motif.set.name, "/"))
  enrich = enrich.test.many.motifs(motif.count.dir, cl)
  save(enrich, file=
    paste0(output.dir, motif.set.name, "/", clustering.name, ".Rdata"))
}

# Looks for motif enrichments for one clustering.
# Args:
#   cl - the clustering
#   name - name to use for the output file
# Side effects: saves motif enrichments to that file.
enrich.test.one.clustering = function(cl, name) {

#   for(motif.set.name in c("Ce_1.02", "Dm_1.02", "Mm_1.02", "Hs_1.02", "bp_1kb_cons0", "meme_1kb_cons0")) {
  for(motif.set.name in c("Ce_1.02", "Dm_1.02", "Mm_1.02", "Hs_1.02")) {
    hyperg.test.1(cl, motif.set.name, name)
  }
}


if (FALSE) {
  # go through the directory of clusterings
  for(clustering.name in c("hier.300.clusters")) {
    # list.files("git/cluster/hierarchical")) {
    cat(clustering.name, "\n")

    cl = as.matrix(read.tsv(
      paste0("git/cluster/hierarchical/", clustering.name, "/clusters.tsv")))[,2]
    cl = gsub(" ", "", cl)

    enrich.test.one.clustering(cl, clustering.name)
  }


}

# Runs on one directory.
enrich.test.file = function(f) {
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

}

if (TRUE) {
  enrich.test.file("hy.tsv")
  enrich.test.file("spencer.tsv")
  enrich.test.file("spencerEmbryonic.tsv") 
#  enrich.test.file("shuffled.rows.facs.tsv")
#  enrich.test.clustering()
}



