# PCA of some samples.

r = as.matrix(read.table("R/unmix/sort_paper/seq/quant/readsPerMillion.tsv",
  header=TRUE, row.names=1, as.is=TRUE))

lr = log2(1+r)

cnd.1.all = lr[,c("cnd.1p8_19", "cnd.1_12_14", "cnd.1p1_4")]

pca.cnd.1 = prcomp(cnd.1.all)

all.reps = lr[,c("cnd.1p8_19", "cnd.1_12_14", "cnd.1p1_4",
  "pha.4p9_1", "pha.4p12_9")]

pca.all.reps = prcomp(all.reps)

pdf("git/unmix/seq/FACS/pca.pdf",
  title="PCA of cnd-1 and pha-4 samples",
  width=7.5, height=7.5)
smoothScatter(pca.all.reps$x[,2], pca.all.reps$x[,3],
  main="PCA of cnd-1 and pha-4 positives",
  xlab="PC2 of genes", ylab="PC3 of genes")

genes.to.show = c(names(sort(pca.all.reps$x[,2])[1:15]),
  names(sort(pca.all.reps$x[,2], decreasing=TRUE)[1:15]),
  names(sort(pca.all.reps$x[,3])[1:15]),
  names(sort(pca.all.reps$x[,3], decreasing=TRUE)[1:15]))
genes.to.show = grep("-", genes.to.show, value=TRUE)
genes.to.show = grep("(pha-4|cnd-1)", genes.to.show, invert=TRUE, value=TRUE)
text(pca.all.reps$x[genes.to.show,2], pca.all.reps$x[genes.to.show,3],
  labels=genes.to.show, cex=0.5)

genes.to.show = c("pha-4", "cnd-1")
text(pca.all.reps$x[genes.to.show,2], pca.all.reps$x[genes.to.show,3],
  labels=genes.to.show, cex=0.9)

plot(pca.all.reps$rotation[,2], pca.all.reps$rotation[,3],
  main="PCA of cnd-1 and pha-4 positives",
  xlab="PC2 rotation of samples", ylab="PC3 rotation of samples",
  type="p", pch=20, xlim=c(-0.8, 0.8), ylim=c(-0.8, 0.8))
text(pca.all.reps$rotation[,2], pca.all.reps$rotation[,3],
  labels=rownames(pca.all.reps$rotation), pos=4)

# comparison with purity numbers from Travis
purity = c(0.94, 0.82, 0.85, 0.96, 0.88)

plot(pca.all.reps$rotation[,2], 100*purity,
  main="Comparison of PC2 with purity measurements",
  xlab="PC2 rotation of sample",
  ylab="FACS purity (%)",
  xlim=c(-0.8, 0.8), ylim=c(80,100), pch=20)
text(pca.all.reps$rotation[,2], 100*purity,
  labels=rownames(pca.all.reps$rotation), pos=4)

plot(pca.all.reps$rotation[,3], 100*purity,
  main="Comparison of PC3 with purity measurements",
  xlab="PC3 rotation of sample",
  ylab="FACS purity (%)",
  xlim=c(-0.8, 0.8), ylim=c(80,100), pch=20)
text(pca.all.reps$rotation[,3], 100*purity,
  labels=rownames(pca.all.reps$rotation), pos=4)

dev.off()

