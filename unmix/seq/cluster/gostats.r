# Attempt at using GOstats.

library("GOstats")
library("org.Ce.eg.db")

# for background genes
source("git/sort_paper/plot/numEnrichedInFractions.r")

# useful utility
gene.to.eg.id = function(genes) {
  mappedLkeys(subset(org.Ce.egALIAS2EG,
    Rkeys = genes, drop.invalid.keys=TRUE))
}

# determine "background" gene symbols
background.genes = names(max.expr)[ max.expr > 1 ]
background.eg.ids = gene.to.eg.id( background.genes )

# Gets GO enrichment results for a set of genes.
# Args:
#   genes - a set of gene symbols (ideally either
#     names like "ceh-6" or clone IDs like "
go.hyperg = function(genes) {
  r = NULL

  for(ont in c("BP", "CC", "MF")) {

    params = new("GOHyperGParams",
      geneIds=gene.to.eg.id(genes),
      universeGeneIds=background.eg.ids,
      annotation="celegans",
      ontology=ont,
      pvalueCutoff=0.05,
      conditional=FALSE,
      testDirection="over")

    ht = summary(hyperGTest(params))
    colnames(ht)[[1]] = "GOID"
    r = rbind(r, cbind(ontology=ont, ht))
  }

  r
}

foo = go.hyperg(sample(background.genes, 100))

