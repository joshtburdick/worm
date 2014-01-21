# Compares FPKM and RPM.

rpm = read.table("git/unmix/seq/quant/readsPerMillion/WS220_20140111/050912.tsv")
rpm.as = read.table("git/unmix/seq/quant/readsPerMillion/WS220_20140111/050912_as.tsv")

gene.fpkm = read.table(gzfile("git/unmix/seq/quant/cufflinks/cufflinks_summary_20140121/Murray_050912_gene_fpkm.tsv.gz"))
isoform.fpkm = read.table(gzfile("git/unmix/seq/quant/cufflinks/cufflinks_summary_20140121/Murray_050912_isoform_fpkm.tsv.gz"))

g = intersect(rownames(rpm), rownames(gene.fpkm))
rpm = rpm[g,]
rpm.as = rpm.as[g,]
gene.fpkm = gene.fpkm[g,]


# one way of computing enrichment
enrich = function(pos, neg, eps=0.001) {
  r = log2(eps + pos) - log2(eps + neg)
  r
}

rpm.enrich = enrich(rpm[,"cehm36p"], rpm[,"cehm36m"])
gene.fpkm.enrich = enrich(gene.fpkm[,"cehm36p"], gene.fpkm[,"cehm36m"])

names(rpm.enrich) = rownames(rpm)
names(gene.fpkm.enrich) = rownames(gene.fpkm)

png("git/unmix/seq/quant/cufflinks/comparisonWithRPM.png",
  width=1000, height=1000)
plot(rpm.enrich, gene.fpkm.enrich, type="p", pch=20, col="#00000080")

dev.off()


a = rpm.enrich - gene.fpkm.enrich
g1 = sort(a, decreasing=TRUE)[1:50]
r1 = cbind(rpm[names(g1),c("cehm36p","cehm36m")],
  gene.fpkm[names(g1),c("cehm36p","cehm36m")])




