# Utilities for converting names, such as gene names.

# WormBase annotation
gene.ids = read.table(gzfile("data/wormbase/c_elegans.PRJNA13758.WS240.xrefs.txt.gz"), sep="\t", as.is=TRUE)[,c(1,2,3,4)]
colnames(gene.ids) = c("gene", "wb.gene", "gene.name", "transcript")
# omit "unknown gene name" cases
gene.ids = gene.ids[ gene.ids$gene.name != "." , ]
gene.id.map = unique(rbind(
  data.frame(id = gene.ids$gene, gene.name=gene.ids$gene.name,
    stringsAsFactors=FALSE),
  data.frame(id = gene.ids$transcript, gene.name=gene.ids$gene.name,
    stringsAsFactors=FALSE)))

# Renames a vector of gene names.
rename.gene.name.vector = function(r) {
  i = r %in% gene.id.map[,"id"]
  r[ i ] =
    gene.id.map[ match(r[i], gene.id.map[,"id"]), "gene.name" ]
  r
}

# Renames the rownames of a table.
# This tries to convert all gene identifiers to gene names,
# where known.
# FIXME this should probably call rename.gene.name.vector().
rename.gene.names = function(a) {
  r = rownames(a)
  i = intersect(r, gene.id.map[,"id"])
  r[ match(i, r) ] =
    gene.id.map[ match(i, gene.id.map[,"id"]), "gene.name" ]

  a = a[ !duplicated(r) , ]
  r = r[ !duplicated(r) ]
  rownames(a) = r
  a
}

