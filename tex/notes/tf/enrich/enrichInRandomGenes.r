# Constructs some random sets of genes, and
# checks for motif enrichment.

source("git/utils.r")
source("git/tf/motif/enrichment/motifHyperg.r")
source("git/sort_paper/FACS/enrichedInFraction.r")

r = read.tsv("git/sort_paper/tf/motif/hyperg/table/facs_0.892.tsv.gz")
ed1 = get.enriched.and.depleted(r.sort.only.averaged, 0.892)

# XXX hokey way of getting this
orig.motif.list = {
  load("git/sort_paper/tf/allResults/motif/hier.300.clusters.Rdata")
  dimnames(enrich)[[1]]
}

random.sets = list()

g = rep(FALSE, nrow(r.sort.only.averaged))
names(g) = rownames(r.sort.only.averaged)

for(n in c(1:15) * 200) {
  a = paste0("random.", n)
  g1 = g
  g1[ sample(names(g1), n) ] = TRUE
  random.sets[[a]] = g1
}

random.enrich = enrich.test.gene.sets.many.motifs(
 "git/tf/motif/count/upstreamMotifCount/5kb/",
 random.sets, orig.motif.list)   

save(random.sets, random.enrich,
  file="git/tex/notes/tf/enrich/enrichInRandomGenes.Rdata")

