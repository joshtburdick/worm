# Finds genes enriched in particular lineages.

source("git/sort_paper/unmix/pseudoinverseEnrichment.r")

# XXX this is somewhat arbitrary
cutoff = 2

x = x.pseudoinverse.23.cell

r = NULL
for(lineage in colnames(x)) {
  lineage = colnames(x)

  x1 = x[,lineage]
  up = sort(x1[ x1 >= cutoff ], decreasing=TRUE)
  down = sort(x1[ x1 <= - cutoff ], decreasing=FALSE)

  if (length(up) > 0)
    r = rbind(r,
      data.frame(set = paste0(lineage, "_enriched"),
        gene = names(up), enrich=up, stringsAsFactors=FALSE))
  if (length(down) > 0)
    r = rbind(r,
      data.frame(set = paste0(lineage, "_depleted"),
        gene = names(down), enrich=down, stringsAsFactors=FALSE)) 
}
rownames(r) = NULL
write.tsv(r, "git/sort_paper/unmix/lineageEnriched.tsv")

