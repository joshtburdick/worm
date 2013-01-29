# Clustering including p-values.

library("pvclust")

r = as.matrix(read.table("git/unmix/seq/quant/readsPerMillion.tsv",
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

colnames(r) = c("ceh-26", "ceh-27", "ceh-36 (-)", "ceh-36 (+)",
  "ceh-6", "cnd-1 (12/14)", "cnd-1 (-)", "cnd-1 (1/4)",
  "cnd-1 (8/19)", "cnd-1 singlets", "cnd-1 ungated", "F21D5.9",
  "hlh-16", "irx-1", "mir-57", "mls-2",
  "pal-1", "pha-4 (-)", "pha-4 (12/9)", "pha-4 (9/1)",
  "pha-4 singlets", "pha-4 ungated", "ttx-3", "unc-130")

r.small = r[1:300,]


p = pvclust(log(1+r), method.dist="uncentered")

pdf("git/unmix/seq/cluster/pvclust.pdf",
  title="clustering of samples", width=10, height=7)
plot(p, cex.pv = 0.8, main=NULL, col.pv = c("red", "white", "white"))
dev.off()

