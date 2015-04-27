# Computes motif enrichment for a few genes.

source("git/utils.r")
source("git/tf/motif/enrichment/motifHyperg.r")
source("git/tf/motif/enrichment/motifHypergSummarize.r")
source("git/tf/motif/motifName.r")

# genes to include
g = read.table("git/tf/pop1/anterior_genes.tsv", as.is=TRUE)[,1]

# where motif counts are
motif.count.base = "git/tf/motif/count/upstreamMotifCount/"

# originally clustered motifs. XXX hokey way of getting this
orig.motif.list = {
  load("git/sort_paper/tf/allResults/motif/hier.300.clusters.Rdata")
  dimnames(enrich)[[1]]
}

orthologs.by.motif.1 = sapply(orthologs.by.motif,
  function(a) paste(a, collapse=" "))


# compute vector for clustering
load(paste0(motif.count.base, "/hughes_20141202/M0103_1.01.Rdata"))

all.genes = dimnames(motif.count)[[1]]
cl = rep(FALSE, length(all.genes))
names(cl) = all.genes
cl[g] = TRUE
gene.sets = list(anterior = cl)

if (FALSE) {
  enrich.hughes = enrich.test.gene.sets.many.motifs(
    paste0(motif.count.base, "/hughes_20141202/"), gene.sets)
  enrich.5kb = enrich.test.gene.sets.many.motifs(
    paste0(motif.count.base, "/5kb/"), gene.sets, motifs=orig.motif.list)
  save(enrich.hughes, enrich.5kb, file="git/tf/pop1/motif_enrich.Rdata")
}

load("git/tf/pop1/motif_enrich.Rdata")

r.hughes = enrich.to.table(enrich.hughes)
r.facs = enrich.to.table(enrich.5kb)

r = rbind(r.hughes, r.facs)
r = r[ order(r$p.corr) , ]
r$orthologs = orthologs.by.motif.1[ r$motif ]
r[ is.na(r$orthologs), "orthologs" ] = ""
write.tsv(r, "git/tf/pop1/anterior_motifs.tsv")

