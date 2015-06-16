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

# annotation of which tissues are in which cell, and coloring
# of them
tissues.per.cell = read.table("data/worm/TissuesPerCell.tsv",
  sep="\t", quote="", header=TRUE, row.names=1, as.is=TRUE)
tissues.per.cell.1 = sapply(lin.node.names,
  function(cell) {
    leaves = names(which(leaf.lineage.matrix[ cell , ] == 1))
    tissues = tissues.per.cell[ leaves, "Tissue" ]
    if (mean(tissues == tissues[1]) == 1)
      tissues[1]
    else
      ""
  })
# various tweaks and simplifications
tissues.per.cell.1[ tissues.per.cell.1 == "Arcade" ] = "Pharynx"
tissues.per.cell.1[ tissues.per.cell.1 == "Nervous" ] = "Neuron"

tissue.to.color = NULL
tissue.to.color["Death"] = "pink"
tissue.to.color["Epidermis"] = "lightblue"
tissue.to.color["Glia"] = "orange"
tissue.to.color["Intestine"] = "darkgreen"
tissue.to.color["Muscle"] = "magenta"
tissue.to.color["Neuron"] = "blue"
tissue.to.color["Pharynx"] = "green"


cell.tissue.colors = tissue.to.color[ tissues.per.cell.1 ]
names(cell.tissue.colors) = lin.node.names
cell.tissue.colors[ is.na(cell.tissue.colors) ] = "grey"

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
    axis(1, at=(j-1)/29, labels=tissue[j],
      col.axis=tissue.colors[j],
      line=0, tick=FALSE, las=2)
  }
  dev.off()
}

plot.one = function(m, hue, rowLabels=NULL) {
  if (is.null(rowLabels))
    rowLabels = rownames(m)

  h = hclust(cor.dist(m))
  m1 = m[h$order,]

  par(mar=c(0,16,0,0.1))
  image(t(m1), useRaster=TRUE,
    xaxt="n", yaxt="n", zlim=c(0,1), col=hsv(hue, 1, 0:128/128))
  n = nrow(m1)-1
  axis(2, at=(0:n)/n, labels=rowLabels, cex.axis=1.39,
    line=0, tick=FALSE, las=2)
}

# Plots these separately.
plot.separately = function() {
  pdf("git/sort_paper/plot/lineage/spencerSortFractionComparison.pdf",
    width=11, height=7)

  layout(matrix(1:3, nrow=3), heights=c(1,1.7,1.3))

  # show the lineage
  par(mar=c(0,16,0,0.1))
  plot(1,1, type="n", xlim=c(24,647), ylim=c(350,0),
    main="", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
  a = rep("#c0c0c0", 1341)
  names(a) = lin.node.names
  par(new=TRUE)
  par(mar=c(0,16,0,0.1))
  plot.segments.per.cell(cell.tissue.colors, main="", root="P0", times=c(0,339),
      lwd=1.5, yaxt="n", int.n.to.label=lin.12.cell, add=TRUE)

  # XXX hack to draw legend outside the plot region
  legend(-100,50, legend=names(tissue.to.color), fill=tissue.to.color,
    bty="n", xpd=NA, cex=1.4)

  plot.one(m.leaf, 0, sapply(rownames(m.leaf), italicize))
  plot.one(spencer.m.leaf, 1/6)
  dev.off()
}

plot.separately()

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


