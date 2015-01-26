# Finds genes enriched in particular lineages.

source("git/sort_paper/unmix/pseudoinverseEnrichment.r")

# XXX this is somewhat arbitrary
cutoff = log2(3)

x = x.pseudoinverse.23.cell

r = NULL
for(j in 1:ncol(x)) {
  lineage = colnames(x)[j]

  x1 = x[,j] - apply(x[,-j], 1, mean)

  up = sort(x1[ x1 >= cutoff ], decreasing=TRUE)
  if (length(up) > 0)
    r = rbind(r,
      data.frame(set = paste0(lineage, "_enriched"),
        gene = names(up), enrich=up, stringsAsFactors=FALSE))

#  down = sort(x1[ x1 <= - cutoff ], decreasing=FALSE)
#  if (length(down) > 0)
#    r = rbind(r,
#      data.frame(set = paste0(lineage, "_depleted"),
#        gene = names(down), enrich=down, stringsAsFactors=FALSE))

}
rownames(r) = NULL
write.tsv(r, "git/sort_paper/unmix/lineageEnriched.tsv")

