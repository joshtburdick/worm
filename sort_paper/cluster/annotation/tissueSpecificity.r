# Plot to estimate of how many clusters are "tissue-specific".

library("VennDiagram")

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

cluster.tissue.spec = apply(abs(cluster.means), 1, mean)
cluster.tissue.spec = cluster.tissue.spec[ names(mean.max.rpm) ]

clusters.to.label = read.tsv("git/sort_paper/cluster/annotation/tissueSpecClustersToLabel.tsv")

pdf("git/sort_paper/cluster/annotation/tissueSpecificity.pdf",
  width=6, height=6)
par(mar=c(5,4,1,1)+0.1)
plot(log2(1+mean.max.rpm), cluster.tissue.spec, pch=20, cex=1,
  col=ifelse(names(mean.max.rpm) %in% anatomy.annotated$cluster,
    "#ff000080", "#00000030"),
  xlab=expression("Mean cluster expression " * (log[2](1 + "RPM"))),
  ylab="Mean absolute enrichment")

cl = rownames(clusters.to.label)
text(log2(1+mean.max.rpm[cl]), cluster.tissue.spec[cl],
    labels = paste0(cl, ": ", clusters.to.label$tissue),
    pos = clusters.to.label$pos, offset = 0.3,
#    adj = c(0,1),
    cex = 0.7)

# showing "tissue-specific" count
rect(4, 0.2, 100, 100, border="#00000040", lwd=2)
dev.off()

# for reference, plot this with all numbers labeled
if (TRUE) {
pdf("git/sort_paper/cluster/annotation/tissueSpecificitySupplemental.pdf",
  width=10, height=10)
a = names(mean.max.rpm) %in% anatomy.annotated$cluster
names(a) = names(mean.max.rpm)
plot(log2(1+mean.max.rpm), cluster.tissue.spec, pch=20, cex=0.5,
  col=ifelse(a, "#ff000080", "#00000030"),
  xlab=expression("Mean cluster expression " * (log[2](1 + "RPM"))),
  ylab="Mean absolute enrichment")
rect(4, 0.2, 100, 100, border="#00000040", lwd=2)
text(log2(1+mean.max.rpm), cluster.tissue.spec,
    labels = names(mean.max.rpm), pos=1,
    col=ifelse(a, "#ff0000a0", "#000000a0"),
    cex = 0.5)
}
dev.off()

# number at various cutoffs
expressed = log2(1 + mean.max.rpm) >= 4
tissue.spec = cluster.tissue.spec >= 0.2
expr.tissue.spec = expressed & tissue.spec

cat("num with tissue spec. >= 0.2 and expr > 4 =",
  sum(expr.tissue.spec), "\n")

cat("proportion of annotated TS above cutoff = ",
  sum(a & expr.tissue.spec), " / ", sum(a), " (",
  round(100 * sum(a & expr.tissue.spec) / sum(a), 2), "%)\n")

cat("number of those annotated =",
  sum(expr.tissue.spec & a[names(expr.tissue.spec)]),
 "\n")

# the Venn diagram
pdf("git/sort_paper/cluster/annotation/vennDiagram.pdf",
  width=6, height=5)

# choose colors
set.seed(43)
h = runif(2)
a1 = hsv(h, 1, 0.5)
a2 = hsv(h, 1, 0.5)

draw.pairwise.venn(sum(expressed), sum(tissue.spec), sum(expr.tissue.spec),
  category = c("Expressed\n\n\n", "Tissue-specific\n\n\n"),
  margin=0.5, lwd=0.5,
  cat.fontfamily = rep("sans", 2),
  cat.col = a1, fill = a2,
  fontfamily = rep("sans", 3))

dev.off()

# table of cluster stats
cluster.stats = data.frame(mean.max.rpm, cluster.tissue.spec, expr.tissue.spec)
rownames(cluster.stats) = names(mean.max.rpm)
write.tsv(cluster.stats,
  "git/sort_paper/cluster/annotation/clusterTissueSpecificity.tsv")

