# Attempt at using GOstats.

library("GOstats")
library("org.Ce.eg.db")

# useful utility
gene.to.eg.id = function(genes) {
  mappedLkeys(subset(org.Ce.egALIAS2EG,
    Rkeys = genes, drop.invalid.keys=TRUE))
}

# Convert parent-link databases to lists, as
# these seem to be faster.
go.parent.list = list(
  BP = as.list(GOBPPARENTS),
  CC = as.list(GOCCPARENTS),
  MF = as.list(GOMFPARENTS))

# Utility to compute depth of a GO term.
# Args:
#   go.ontology - which GO ontology ("BP", "CC", or "MF")
#   go.id - a GO ID
# Returns: depth of term with that GO ID
go.depth = function(go.ontology) {
  go.parent.db = go.parent.list[[ go.ontology ]]

  function(go.id) {
    i = 0
    s = c(go.id)

    # simple version of breadth-first search
    while (!("all" %in% s)) {
      s1 = c(as.vector(sapply(s, function(x) as.vector(go.parent.db[[x]]))),
        recursive=TRUE)
      s = unique(c(s, s1))
      i = i + 1
    }

    return(i)
  }
}

# Computes depth for many GO terms.
# Args:
#   go.category - either "BP", "CC", or "MF"
#   go.id - a vector of GO IDs
# Returns: depth of those IDs
go.depth.deprecated = function(go.category, go.id) {
  parent.db = NULL
  if (go.category == "BP") parent.db = GOBPPARENTS
  if (go.category == "CC") parent.db = GOCCPARENTS
  if (go.category == "MF") parent.db = GOMFPARENTS

  f = function(x) go.depth.1( , x)

  sapply(go.id, f)
}

# Gets GO enrichment results for a set of genes.
# Args:
#   background.genes - the set of background genes
#   genes - a set of gene symbols (ideally either
#     names like "ceh-6" or clone IDs like "T19A6.1")
go.hyperg = function(background.genes, genes) {
  background.eg.ids = gene.to.eg.id( background.genes )

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

    try({
      ht = summary(hyperGTest(params))
      colnames(ht)[[1]] = "GOID"
      print(dim(ht))
      if (nrow(ht) > 0) {
        r = rbind(r, cbind(ontology=ont, ht))
      }
    })
  }

  r$p.corr = p.adjust(r$Pvalue, method="fdr")
  r = r[r$p.corr <= 0.05,]
  r
}



