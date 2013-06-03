# Gets eigengenes from WGCNA analysis.

options(stringsAsFactors = FALSE)

source("git/unmix/seq/cluster/WGCNA/wgcna.r")

load("git/unmix/seq/cluster/WGCNA/wnet.ts/wnet.Rdata")
wnet = wnet.ts    # XXX

# moduleColors = labels2colors(wnet$colors)

module.eigengene =
  moduleEigengenes(t(r), wnet$colors, trapErrors=TRUE)$eigengenes
rownames(module.eigengene) = colnames(r)
module.eigengene = as.matrix(module.eigengene)

# sort this so that clusters resemble ordering of fractions
# s = t(module.eigengene) %*% c(1:64)^2
# s = t(apply(module.eigengene, 2, order)) %*% c(1:64)
f = function(x) {
  o = order(x, decreasing=TRUE)
  o[1] + 0.1 * o[2] + 0.01 * o[3]
}
s = apply(module.eigengene, 2, f)

module.eigengene = module.eigengene[,order(s)]
colnames(module.eigengene) = c(1:ncol(module.eigengene))

# if (TRUE) {
  # table of what the eigengenes are
  write.table(round(module.eigengene, 5),
    sep="\t", col.names=NA,
      file="git/unmix/seq/cluster/WGCNA/eigengene.ts.txt")

  # degree of membership of each gene in each cluster
  kme = as.matrix(signedKME(t(r), module.eigengene))

  # table of which gene is in which cluster
  # FIXME add degree-of-membership?
  gene.cluster = data.frame(
    cluster = order(order(s))[wnet$color+1],
#orig.cluster = wnet$color,    
    row.names = rownames(r))
  gene.cluster$kme =
    round(kme[ cbind(c(1:nrow(kme)), gene.cluster$cluster) ], 3)
  gene.cluster = gene.cluster[
    order(gene.cluster$cluster, -gene.cluster$kme), ]
  write.table(gene.cluster, sep="\t", col.names=NA,
    file="git/unmix/seq/cluster/WGCNA/geneCluster.tsv")

# }

