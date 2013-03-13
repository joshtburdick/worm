# Analyzes groups of genes with GOstats.

library("AnnotationForge")
library("org.Ce.eg.db")
library("GSEABase")
library("GOstats")
library("biomaRt")

f = toTable(org.Ce.egGO)

goframeData = data.frame(f$go_id, f$Evidence, f$gene_id)

goFrame = GOFrame(goframeData,organism="Caenerhabditis elegans")
goAllFrame = GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

# Biomart connections.
wb.gene = useMart("WS220", "wormbase_gene")
ensembl = useMart("ensembl", dataset="celegans_gene_ensembl")

# Generalized "conversion" function.
# Args:
#   mart - a biomaRt object
#   from - attribute to filter on
#   to - attribute to return
# Returns: function which converts names.
mart.convert = function(mart, from, to) function(x) {
  r = getBM(attributes = c(to),
    filters = c(from),   
    values = list(x),
    mart = mart)
  unique(r)
}

# Useful conversion functions.
wb.gene.name.to.gene.id =
  mart.convert(wb.gene, "gene_name", "sequence_name")
gene.id.to.entrez = 
  mart.convert(ensembl, "ensembl_gene_id", "entrezgene")


# Computes enriched GO categories.
# Args:
#   gene - the gene names
#   cluster - the cluster each gene is in
go.enrichment = function(gene, cluster) {

cat("getting IDs for all genes...")
  universe = gene.id.to.entrez(wb.gene.name.to.gene.id(gene))
cat("done\n")

  go = list(BP=NULL, CC=NULL, MF=NULL)

  for(cl in sort(unique(cluster))) {
    g = gene[ cluster == cl ]
    gene.list = gene.id.to.entrez(wb.gene.name.to.gene.id(g))

    for(ont in c("BP", "CC", "MF")) {
cat(cl, ont, "  ")
      params = GSEAGOHyperGParams(name="HyperG params",
        geneSetCollection=gsc,
        geneIds = gene.list,
        universeGeneIds = universe,
        ontology = ont,
        pvalueCutoff = 0.05,
        conditional = FALSE,
        testDirection = "over")
      r = hyperGTest(params)
      s = summary(r)
      if (nrow(s) > 0)
        go[[ont]] = rbind(go[[ont]],
          cbind(cluster=as.character(cl), s))
    }
  }

  go
}


# gene.id.to.entrez(c("171590", "171591", "171592"))
# gene.id.to.entrez(c("T24D1.1"))


clustering = read.table(
  "git/unmix/seq/cluster/WGCNA/geneCluster.tsv",
  sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)

go = go.enrichment(rownames(clustering), clustering$cluster)

save(go, file="git/unmix/seq/cluster/WGCNA/gostats.Rdata")
write.table(go[["BP"]], file="gostats/BP.tsv", sep="\t", row.names=TRUE, col.names=NA)
write.table(go[["CC"]], file="gostats/CC.tsv", sep="\t", row.names=TRUE, col.names=NA)
write.table(go[["MF"]], file="gostats/MF.tsv", sep="\t", row.names=TRUE, col.names=NA)

