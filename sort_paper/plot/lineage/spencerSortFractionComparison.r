# Clusters our sort fractions with the Spencer
# sort fractions.

source("git/utils.r")
source("git/plot/utils.r")
source("git/plot/plot_expr.r")

# the sort matrix (unnormalized)
source("git/sort_paper/unmix/sortMatrix.r")

# just show the leaf nodes
m.leaf = m.unnormalized %*% cell.to.leaf.matrix

# XXX
colnames(m.leaf)[colnames(m.leaf)=="NA."] = "P0"

# omitting double-negatives, for consistency with other rows
m.leaf = m.leaf[ rownames(m.leaf) != "ceh-6 (-) hlh-16 (-)" , ]

spencer.m = as.matrix(read.tsv(
  "git/sort_paper/validation/spencerSortMatrix.tsv"))
spencer.m.leaf = spencer.m %*% cell.to.leaf.matrix

m1 = rbind(m.leaf, spencer.m.leaf)
m1 = 2 * m1
m1[ m1 > 1 ] = 1

# definition of correlation, which skips missing values
cor.dist = function(a) {
  a = as.dist( 1 - cor(t(as.matrix(a)), use="na.or.complete") )
#  a[ is.na(a) ] = 2
  a
}


plot.it.orig = function() {
  h = hclust(cor.dist(m1))

  pdf("git/sort_paper/plot/lineage/spencerSortFractionComparison.pdf",
    width=7, height=10)

  layout(t(t(c(1,2))), heights=c(1,3))

  par(mar=c(0,0,1,0))
  plot(h, labels=FALSE, axes=FALSE,
    sub="", main="", xlab="", ylab="", hang=-2)

  par(mar=c(11,0.7,0,0.7))
  image(m1[h$order,], useRaster=TRUE,
    xaxt="n", yaxt="n", zlim=c(0,1), col=hsv(1/6, 1, 0:128/128))

  tissue = h$labels[h$order]
  tissue.colors = ifelse(tissue %in% rownames(m.leaf), "red", "blue")

  # XXX somewhat hacky
  for(j in 1:30) {
    axis(1, at=(j-1)/29, labels=tissue[j], col.axis=tissue.colors[j],
      line=0, tick=FALSE, las=2)
  }
  dev.off()
}

plot.one = function(m, col) {
  h = hclust(cor.dist(m))

  plot(h, labels=FALSE, axes=FALSE,
    sub="", main="", xlab="", ylab="", hang=-2)


}



# Various comparisons of the sizes of these groups.
facs.on.off = m.leaf > 0
spencer.on.off = spencer.m.leaf > 0

mean.overlap = function(m) {
  a = m %*% t(m)
  diag(a) = NA
  a = as.vector(a)
  mean(as.vector(a), na.rm=TRUE)
}


cat(paste0("FACS num cells    mean = ",
  mean(apply(facs.on.off,1,sum)),
  "    median = ", median(apply(facs.on.off,1,sum)), "\n"))
cat(paste0("Spencer num cells    mean = ",
  mean(apply(spencer.on.off,1,sum)),
  "    median = ", median(apply(spencer.on.off,1,sum)), "\n"))

cat(paste0("FACS num. fractions    mean = ",
  mean(apply(facs.on.off,2,sum)),
  "    median = ", median(apply(facs.on.off,2,sum)), "\n"))
cat(paste0("Spencer num. fractions   mean = ",
  mean(apply(spencer.on.off,2,sum)),
  "    median = ", median(apply(spencer.on.off,2,sum)), "\n"))

cat(paste0("FACS mean overlap = ", mean.overlap(facs.on.off), "\n"))
cat(paste0("Spencer mean overlap = ", mean.overlap(spencer.on.off), "\n"))


