# Utilities for converting names, such as gene names.

gene.ids = read.csv(gzfile("data/wormbase/geneIDs.WS224.csv.gz"),
  as.is=TRUE, header=FALSE)
colnames(gene.ids) = c("wb.gene", "gene", "transcript")
gene.ids = gene.ids[ gene.ids[,"gene"] != "" , ]

# Renames the row names of a matrix (or data frame).
rename.gene.names = function(a) {
  r = rownames(a)
  i = intersect(r, gene.ids[,"transcript"])
  r[ match(i, r) ] =
    gene.ids[ match(i, gene.ids[,"transcript"]), "gene" ]

  a = a[ !duplicated(r) , ]
  r = r[ !duplicated(r) ]
  rownames(a) = r
  a
}

