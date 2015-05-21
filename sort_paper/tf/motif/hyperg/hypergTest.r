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
    "/home/jburdick/gcb/git/sort_paper/tf/motif/hughes/upstreamCount/"
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

#  for(motif.set.name in c("Mm_1.02", "Hs_1.02", "bp_1kb_cons0", "meme_1kb_cons0")) {
  for(motif.set.name in c("Ce_1.02")) {
    hyperg.test.1(cl, motif.set.name, name)
  }

# everything else is deprecated
  # e.g., this is not currently used
  output.dir = "git/sort_paper/tf/motif/hyperg/allResults_check/"

if (FALSE) {
  system(paste0("mkdir -p ", output.dir, "chip/"))
  enrich = enrich.test.many.motifs(
    "git/tf/chip/count/upstreamChipCount/", cl)
  save(enrich, file=
    paste0(output.dir, "chip/", name, ".Rdata"))
}
if (FALSE) {
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
  for(clustering.name in c("hier.300.clusters")) {    # list.files("git/cluster/hierarchical")) {
    cat(clustering.name, "\n")

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



