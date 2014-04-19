# Plots the expression patterns of the markers
# used for sorting, and the cluster centers.

source("git/utils.r")

# the sort matrix (unnormalized)
source("git/sort_paper/unmix/sortMatrix.r")

m.leaf = m.unnormalized %*% cell.to.leaf.matrix

# one of the clusterings
clustering = {
  cl1 = read.table("git/cluster/hierarchical/hier.50.clusters/clusters.tsv",
    sep="\t", header=TRUE, as.is=TRUE)
  cl = cl1[,3]
  names(cl) = cl1[,2]
  cl
}

# the enrichments
r = as.matrix(read.tsv("git/cluster/readRatios.tsv"))
r = r[ names(cl) , ]
r = r[ names(clustering) , ]

# cluster centers
cluster.mean = t(apply(r, 2, function(x)
  c(by(x, clustering, function(y) mean(y, na.rm=TRUE)))))

cluster.mean = rbind( cluster.mean[ c(23:1, 38:24), ])


#pdf("git/sort_paper/plot/heatmap/sortMarkersAndClusters.pdf",
#  width=12, height=8)
png("git/sort_paper/plot/heatmap/sortMarkersAndClusters.png",
  width=1200, height=600)

layout(matrix(1:4, nrow=2),
  widths=c(2,3), heights=c(1,3))

h = hclust(dist(t(cluster.mean)))

plot.new()
marker.on.off = rbind(m.leaf[c(1:12,11,11,11,4,4,4),],
  m.blank)


image(t(marker.on.off),
  xaxt="n", yaxt="n", zlim=c(0,1),
  col=hsv(0, 1, 0:128/128))
plclust(h, labels=FALSE, axes=FALSE)

image(t(cluster.mean[,h$order]))

dev.off()

