# Attempt at using GOstats. Not currently working.


library("GOstats")
library("org.Ce.eg.db")

library("biomaRt")   # this is deprecated

ens.ce = useMart("ensembl", dataset="celegans_gene_ensembl")

# to see available maps, do
#    ls("package:celegans.db")

# convert symbol names to Entrez gene IDs
symbol.to.entrez = function(genes) {
  a = unlist(mget(genes, revmap(celegansSYMBOL), ifnotfound=NA))
  a = a[ !is.na(a) ]
  s = unlist(mget(a, celegansENSEMBL, ifnotfound=NA))
  s = s[ !is.na(s) ]
  s
}

# print( symbol.to.entrez( c("ceh-26", "scl-27", "sre-20", "foo") ))
gene.ids = c("ceh-6", "hlh-1", "tbx-35", "cnd-1")

all.genes = Lkeys(org.Ce.egGO)

gene.ids = sample(all.genes, 500)

# parameters
if (TRUE) {
params = new("GOHyperGParams",
  geneIds=gene.ids,
#  universeGeneIds= symbol.to.entrez(c("lin-12", "lin-14")),
  annotation="celegans",
  ontology="BP",
  pvalueCutoff=0.05,
  conditional=FALSE,
  testDirection="over")
}

r = hyperGTest(params)

