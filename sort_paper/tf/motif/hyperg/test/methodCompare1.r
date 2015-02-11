# Compares the clustering-based and set-based methods of
# looking for motif enrichments.

source("git/tf/motif/enrichment/motifHyperg.r")

# for comparison's sake
load("git/sort_paper/tf/motif/hyperg/allResults/5kb/hier.300.clusters.Rdata")

# XXX hokey way of getting this
orig.motif.list = {
  load("git/sort_paper/tf/allResults/motif/hier.300.clusters.Rdata")
  dimnames(enrich)[[1]]
}

# a clustering
clustering = "hier.300.clusters"
cl1 = read.tsv(paste0("git/cluster/hierarchical/", clustering, "/clusters.tsv"))
cl = cl1[,2]
names(cl) = rownames(cl1)

# Converts from a clustering to a list of boolean vectors;
# possibly I wrote this elsewhere.
cl.to.vectors = function(cl) {
  r = list()
  for (a in as.character(sort(unique(cl)))) {
    r[[ a ]] = (cl == a)
  }
  r
}

gene.sets = cl.to.vectors(cl)

enrich1 = enrich.test.gene.sets.many.motifs(
  "git/tf/motif/count/upstreamMotifCount/5kb/",
  gene.sets, orig.motif.list)

# compare all the stats
print(apply(enrich - enrich1, c(3), range))

save(enrich1,
  file="git/sort_paper/tf/motif/hyperg/test/methodCompare1.Rdata")

