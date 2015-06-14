# Plot to estimate of how many clusters are "tissue-specific".

source("git/utils.r")

cl = as.matrix(read.tsv(
  paste0("git/cluster/hierarchical/hier.300.clusters/clusters.tsv")))[,2]
cl = gsub(" ", "", cl)

cluster.sizes = table(cl)

cat(paste0("cluster sizes ranged from ", min(cluster.sizes),
  " to ", max(cluster.sizes), " genes\n"))

# clusters with some anatomy annotation enriched
anatomy.annotated = read.tsv(
  "git/sort_paper/enrichment/clustersWithAnatomyAnnotation.tsv")

# for comparing cluster centers and TF expression profiles
cluster.means.1 = read.tsv(
  paste0("git/sort_paper/cluster/centroids/hier.300.clusters_means.tsv"))
cluster.means = as.matrix(cluster.means.1[,1:23])

# reads per million (for computing average of maximum)
rpm = read.tsv("git/cluster/readsPerMillion.tsv")
rpm = rpm[ names(cl), !grepl("^(HS|RNAi)", colnames(rpm)) ]
max.rpm = apply(rpm, 1, max)
mean.max.rpm = tapply(max.rpm, cl, mean)
mean.max.rpm = mean.max.rpm[ order(as.numeric(names(mean.max.rpm))) ]

cluster.ratio.sd = apply(cluster.means, 1, sd)
cluster.ratio.sd = cluster.ratio.var[ names(mean.max.rpm) ]

clusters.to.label = read.tsv("git/sort_paper/cluster/annotation/tissueSpecClustersToLabel.tsv")

pdf("git/sort_paper/cluster/annotation/tissueSpecificity.pdf",
  width=6, height=5)

plot(log2(1+mean.max.rpm), cluster.ratio.sd, pch=20, cex=1,
  col=ifelse(names(mean.max.rpm) %in% anatomy.annotated$cluster,
    "#ff000080", "#00000030"),
  xlab="Mean cluster expression (log2(1 + RPM))",
  ylab="Standard deviation of cluster enrichments")

cl = rownames(clusters.to.label)
text(log2(1+mean.max.rpm[cl]), cluster.ratio.sd[cl],
    labels = paste0(cl, ": ", clusters.to.label$tissue),
    pos = clusters.to.label$pos, offset = 0.3,
#    adj = c(0,1),
    cex = 0.7)

# showing "tissue-specific" count
rect(4, 0.1, 100, 100, border="#00000040", lwd=2)

# for debugging
if (FALSE) {
text(log2(1+mean.max.rpm), cluster.ratio.sd,
    labels = names(mean.max.rpm), pos=1,
    cex = 0.5)
}
dev.off()


# number at various cutoffs
cat("num with sd >= 0.1 and expr > 4 =",
  sum(((log2(1 + mean.max.rpm) >= 4) & (cluster.ratio.sd >= 0.1))),
  "\n")




