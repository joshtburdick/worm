# Plots the expression patterns of the markers
# used for sorting, and the cluster centers.

source("git/utils.r")
source("git/plot/label_panel.r")
source("git/plot/utils.r")
source("git/plot/plot_expr.r")

# the sort matrix (unnormalized)
source("git/sort_paper/unmix/sortMatrix.r")

m.leaf = m.unnormalized %*% cell.to.leaf.matrix

# one of the clusterings
clustering = {
  cl1 = read.table("git/cluster/hierarchical/hier.300.clusters/clusters.tsv",
    sep="\t", header=TRUE, as.is=TRUE)
  cl = cl1[,3]
  names(cl) = cl1[,2]
  cl
}

# the enrichments
r = as.matrix(read.tsv("git/cluster/readRatios.tsv"))
r = r[ names(cl) , ]
r = r[ names(clustering) , ]

# average the enrichments for cnd-1 and pha-4
r = cbind(
  "cnd-1"=apply(r[,c("cnd-1 8/19","cnd-1 12/14","cnd-1 1/4")], 1, mean),
  "pha-4"=apply(r[,c("pha-4 5/9","pha-4 9/1","pha-4 12/9")], 1, mean),
  r[ , c(7:38) ])

# cluster centers
cluster.mean = t(apply(r, 2, function(x)
  c(by(x, clustering, function(y) mean(y, na.rm=TRUE)))))

# XXX somewhat hacky row sorting
fraction.names = c("ceh-36", sort(setdiff(colnames(r)[1:14], "ceh-36")))

cluster.mean = rbind(
  cluster.mean[ c(19:15) , ],
  cluster.mean[ rev(fraction.names) , ],
  cluster.mean[ c(34:20) , ])


# where markers are on and off
marker.on.off = matrix(NA, nrow=nrow(cluster.mean), ncol=671)
rownames(marker.on.off) = rownames(cluster.mean)
colnames(marker.on.off) = colnames(m.leaf)

# only include samples with image and expression data
m1 = intersect(rownames(m.leaf), rownames(cluster.mean))
marker.on.off[ m1 , ] = m.leaf[ m1 , ]
# cluster.mean = cluster.mean[ m1 , ]

# FIXME: use the original clustering
h = hclust(dist(t(cluster.mean[1:20,])))

pdf("git/sort_paper/plot/heatmap/sortMarkersAndClusters.pdf",
  width=16, height=8)
# png("git/sort_paper/plot/heatmap/sortMarkersAndClusters.png",
#   width=1200, height=600)
layout(matrix(c(1,2,3,4,4,5), nrow=2, byrow=TRUE),
  widths=c(1,1,3), heights=c(1,4))

blue.yellow.colors =
#  c(hsv(2/3, 1, 128:0/128), hsv(1/6, 1, 1:128/128))
# temporarily trying to match TreeView colors
  c(hsv(194/360, 1, 128:0/128), hsv(60/360, 1, 1:128/128))

# color bar for markers
par(mar=c(4,2,7,6)+0.1)
image(as.matrix(0:128/128), col=hsv(0, 1, 0:128/128),
  xaxt="n", yaxt="n", useRaster=TRUE, main="Marker expression")
axis(1, label=c("off","on"), at=c(0,1))

# color bar for expression
par(mar=c(4,5,7,3)+0.1)
image(as.matrix(-128:128/128), col=blue.yellow.colors,
  xaxt="n", yaxt="n", useRaster=TRUE,
  main="Enrichment measured by RNA-seq")
axis(1, label=c("depleted", "no difference", "enriched"), at=c(0,0.5,1))

# clustering dendrogram
par(mar=c(0,4,0.1,0))
plclust(h, labels=FALSE, axes=FALSE,
  sub="", xlab="", ylab="", hang=-1)

# ??? brighten image a bit?
# marker.on.off[ marker.on.off >= 0.5 ] = 1

# show where markers are on and off
par(mar=c(2.1,2.1,0,2.1))
image(t(marker.on.off),
  useRaster=TRUE,
  xaxt="n", yaxt="n", zlim=c(0,1),
  col=hsv(0, 1, 0:128/128))

# draw the ceh-36 tree above this
# XXX this is really hacky
par(new=TRUE)
plot(1,1, type="n", xlim=c(22,649), ylim=c(1250,40),
  xaxt="n", yaxt="n")
par(new=TRUE)
a = rgb(m.unnormalized["ceh-36",], 0, 0)
names(a) = lin.node.names
plot.segments.per.cell(a, main="", root="P0", times=c(0,550),
    lwd=1.5, yaxt="n", add=TRUE)

label.panel("a)", gp=gpar(fontsize=14, col="black"))


# the actual clustered expression
par(mar=c(2.1,6.4,0,2.4))
image(scale.interval.to.unit(t(cluster.mean[,h$order]), c(-1,1)),
  useRaster=TRUE, 
  xaxt="n", yaxt="n", zlim=c(0,1),
  col=blue.yellow.colors)
label.panel("b)", gp=gpar(fontsize=14, col="black"))

# label the samples
sample.labels = rownames(cluster.mean)
sample.labels = sub("t\\.00", "", sample.labels)
sample.labels = sub("t\\.0", "", sample.labels)
sample.labels = sub("t\\.", "", sample.labels)
axis(2, at=0:33/33, labels=sample.labels,
  las=2, tick=FALSE, line=-0.7)

mtext("Time (minutes)", side=2, adj=0.87, line=5.5, cex=1.3)
mtext("FACS sample", side=2, adj=0.3, line=5.5, cex=1.3)

dev.off()

