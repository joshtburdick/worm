# Computes motif enrichment for a few genes.

source("git/utils.r")
source("git/util/arrayUtils.r")
source("git/tf/motif/enrichment/motifHyperg.r")
source("git/tf/motif/enrichment/motifHypergSummarize.r")
# source("git/tf/motif/motifName.r")
source("git/sort_paper/tf/motif/hughes/motifInfo.r")

# genes to include
g.anterior = read.table("git/tf/pop1/anterior_genes.tsv", as.is=TRUE)[,1]
g.posterior = read.table("git/tf/pop1/posterior_genes.tsv", as.is=TRUE)[,1]

# where motif counts are
motif.count.base = "git/sort_paper/tf/motif/upstreamCount/"

orthologs.by.motif.1 = sapply(orthologs.by.motif,
  function(a) paste(a, collapse=" "))

# Converts from a list of genes to a boolean vector.
# (This is just to get the list of all genes for which we have motifs.)
load(paste0(motif.count.base, "/Ce_1.02/M0103_1.02.Rdata"))
genes.to.cluster = function(g) {
  all.genes = dimnames(motif.count)[[1]]
  cl = rep(FALSE, length(all.genes))
  names(cl) = all.genes
  cl[g] = TRUE
  cl
}
gene.sets = list(anterior = genes.to.cluster(g.anterior),
  posterior = genes.to.cluster(g.posterior))

# where to store results
output.dir = "git/tf/pop1/enrichResults/"

# Test for enriched motifs.
if (FALSE) {

  for(motif.set in c("Ce_1.02", "Dm_1.02", "Mm_1.02", "Hs_1.02")) {
    cat(motif.set, "\n")
    system(paste0("mkdir -p ", output.dir, "/", motif.set))
    enrich = enrich.test.gene.sets.many.motifs(
      paste0(motif.count.base, "/", motif.set), gene.sets)
    save(enrich, file=paste0(output.dir, "/", motif.set, "/pop1.Rdata"))
  }
}


# Convert to a table.
if (TRUE) {

  # read in results from each organism
  enrich1 = list()
  for(motif.set in c("Ce_1.02", "Dm_1.02", "Mm_1.02", "Hs_1.02")) {
    load(paste0(output.dir, "/", motif.set, "/pop1.Rdata"))
    enrich1[[ motif.set ]] = enrich
  }

  r = enrich.to.table.many(enrich1)

  # convert to a table, correcting p-values of combined results
  r = r[ order(r$p.corr) , ]

  # add some annotation
  r$motif.name = motif.info[r$motif, "motif.name"]
  r = cbind(r, motif.info[r$motif, c("species", "related.gene")])
  r$orthologs = orthologs.by.motif.1[ r$motif ]

  # for Ce genes, count the actual gene as an "ortholog"
  i = (r$species=="C.elegans")
  r[ i , "orthologs" ] = r[ i, "related.gene" ]

  write.tsv(r, "git/tf/pop1/motif_enrich.tsv")
}


