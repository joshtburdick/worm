# Genes enriched in individual fractions.

output.dir = "git/unmix/seq/cluster/"

r = as.matrix(read.table(
    "git/unmix/seq/cluster/readsNormalized.tsv",
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

le = r

enriched.and.depleted =
  data.frame(enriched = apply(le>=2, 2, sum),
    depleted=apply(le <=-2, 2, sum),
    stringsAsFactors=FALSE)

write.table(enriched.and.depleted, sep="\t",
  file=paste(output.dir, "/enriched.and.depleted_.tsv", sep=""))

# write directory of files of which genes are enriched in each fraction
system(paste("mkdir -p ", output.dir, "/enriched_genes_2/", sep=""))
for(g in colnames(le)) {
  enriched = rownames(le)[ le[,g] >= 2 ]
  enriched = gsub(" ", "_", enriched)
  enriched = gsub("/", "_", enriched)
  write.table(enriched, quote=FALSE, row.names=FALSE, col.names=FALSE,
    file=paste(output.dir, "/enriched_genes/", g, ".tsv", sep=""))
}

