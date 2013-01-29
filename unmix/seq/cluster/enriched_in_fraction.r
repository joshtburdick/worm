# Genes enriched in individual fractions.

output.dir = "git/unmix/seq/cluster/"

r = as.matrix(read.table(gzfile(
    "git/unmix/seq/quant/readsPerMillion_pooled.tsv.gz"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

lr = log2(1 + r)
le = lr - lr[,"all"]

# tack on cases with negatives
genes.with.neg = c("ceh-36", "cnd-1", "pha-4")
le.neg = le[,genes.with.neg] - le[,paste(genes.with.neg, "_minus", sep="")]
colnames(le.neg) = paste(genes.with.neg, "_rel_to_neg", sep="")
le = cbind(le, le.neg)

enriched.and.depleted =
  data.frame(enriched = apply(le>=2, 2, sum),
    depleted=apply(le <=-2, 2, sum),
    stringsAsFactors=FALSE)

write.table(enriched.and.depleted, sep="\t",
  file=paste(output.dir, "/enriched.and.depleted.tsv", sep=""))

# write directory of files of which genes are enriched in each fraction
system(paste("mkdir -p ", output.dir, "/enriched_genes/", sep=""))
for(g in colnames(le)) {
  enriched = rownames(le)[ le[,g] >= 2 ]
  write.table(enriched, quote=FALSE, row.names=FALSE, col.names=FALSE,
    file=paste(output.dir, "/enriched_genes/", g, ".tsv", sep=""))
}





