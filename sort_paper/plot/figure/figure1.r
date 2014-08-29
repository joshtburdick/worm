# Figure 1: showing lineages used for sorting, and the
# corresponding scatterplots.


source("git/unmix/seq/plot/scatterplot.r")

source("git/sort_paper/unmix/sortMatrix.r")
source("git/plot/plot_expr.r")
source("git/plot/label_panel.r")

# font size for panel labels
gp = gpar(fontsize=12, col="black")

lineage.to.plot = "ABplpa"
cells = get.node.names(get.lineage.by.name(lin, lineage.to.plot), 100)
cells = colnames(m.unnormalized)[ colnames(m.unnormalized) %in% cells ]

m1 = m.unnormalized[ c("ceh-27", "hlh-16", "unc-130", "irx-1"), cells ]

# Plots a lineage in a particular color.
plot.lineage.small = function(gene, main, hue) {

#  par(mar=c(4,4,3,3)+0.1)
  par(mar=c(0,0,0,0))
  a = hsv(hue, 1, 0.8 * m.unnormalized[gene,])

  names(a) = lin.node.names

  plot.segments.per.cell(a, main, root="ABplpa", times=c(50,350),
      lwd=6, yaxt="n")
  axis(2, font=1)
#  mtext("Time (minutes)", side=2, line=2, cex=0.7)

}

pdf("git/sort_paper/plot/figure/figure1.pdf",
  width=8, height=8)

layout(matrix(c(1:10), nrow=5),
  widths=c(2,1), heights=c(1,1,1,1,2))

# plot the trees
plot.lineage.small("ceh-27", "", 0)
plot.lineage.small("hlh-16", "", 1/6)
plot.lineage.small("unc-130", "", 2/6)
plot.lineage.small("irx-1", "", 3/6)

# made-up expression
x = matrix(0, nrow=ncol(m1), ncol=5)
x[1:20, 1:2] = 1
x[5:8, 3] = 1
x[34:49, 4:5] = 1

x = x+rnorm(prod(dim(x)), mean=0, sd=0.2)
x[x < 0] = 0
x[x > 1] = 1

par(mar=c(0,0,1,0) + 0.1)
image(x, col=hsv(0,0,128:0/128), useRaster=TRUE,
  xaxt="n", yaxt="n")

# utility to color a slice of the heatmap
color.slice = function(y1, y2, hue) {
  rect(-1, (y1-0.5)/4, 2, (y2-0.5)/4, col=hsv(hue, 1, 1, 0.1), border=NA)
}

color.slice(0, 2, 0.1)
color.slice(2, 3, 0.3)
color.slice(3, 5, 0.7)

plot.new()
plot.new()
plot.new()
plot.new()

x.p = m1 %*% x
x.n = (1-m1) %*% x
x.f = x.p / (x.p + x.n)

par(mar=c(0,1,1,0) + 0.1)
image(x.f, col=hsv(0,0,128:0/128), useRaster=TRUE,
  xaxt="n", yaxt="n")
color.slice(0, 2, 0.1)
color.slice(2, 3, 0.3)
color.slice(3, 5, 0.7)

dev.off()

