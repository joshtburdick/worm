# Genes enriched in individual fractions.

output.dir = "git/unmix/seq/cluster/"

r = as.matrix(read.table(
    "git/unmix/seq/cluster/readsNormalized.tsv",
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

le = r

enrich.cutoff = log2(4)

enriched.and.depleted =
  data.frame(enriched = apply(le>=enrich.cutoff, 2, sum),
    depleted=apply(le <=-enrich.cutoff, 2, sum),
    stringsAsFactors=FALSE)

write.table(enriched.and.depleted, sep="\t",
  file=paste(output.dir, "/enriched.and.depleted.tsv", sep=""),
  row.names=TRUE, col.names=NA)

# write directory of files of which genes are enriched in each fraction
system(paste("mkdir -p ", output.dir, "/enriched_genes/", sep=""))
for(g in colnames(le)) {
  enriched = rownames(le)[ le[,g] >= enrich.cutoff ]
  g = gsub(" ", "_", g)
  g = gsub("/", "_", g)
  write.table(enriched, quote=FALSE, row.names=FALSE, col.names=FALSE,
    file=paste(output.dir, "/enriched_genes/", g, ".tsv", sep=""))
}

