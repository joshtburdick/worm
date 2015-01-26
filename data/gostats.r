# Attempt at using GOstats.

library("GOstats")
library("GO.db")
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
#   universe.genes - the universe of gene names
#   genes - a set of gene symbols (ideally either
#     names like "ceh-6" or clone IDs like "T19A6.1")
#   depth.bound - only include terms with depth between
#     depth.bound[1] and depth.bound[2]
go.hyperg = function(universe.genes, genes, depth.bound) {
  universe.eg.ids = unique(gene.to.eg.id( c(universe.genes, genes)))

  r = NULL
  num.tests = 0

  for(ont in c("BP", "CC", "MF")) {
    geneIds = gene.to.eg.id(genes)

    if (length(gene.to.eg.id) > 0) {
      params = new("GOHyperGParams",
        geneIds=gene.to.eg.id(genes),
        universeGeneIds=universe.eg.ids,
        annotation="celegans",
        ontology=ont,
        pvalueCutoff=0.05,
        conditional=TRUE,   # was FALSE
        testDirection="over")

      try({
        gt = hyperGTest(params)
        ht = summary(gt, categorySize=10)
        colnames(ht)[[1]] = "GOID"
        print(dim(ht))
        num.tests = num.tests + length(slot(slot(gt, "goDag"), "nodes"))
        if (nrow(ht) > 0) {
          r = rbind(r, cbind(ontology=ont, ht,
            term.depth = sapply(ht$GOID, go.depth(ont))))
        }
      })
    }
  }

  r = r[ r$term.depth >= depth.bound[1] & r$term.depth <= depth.bound[2] , ]
  list(r = r, num.tests = num.tests)
}

