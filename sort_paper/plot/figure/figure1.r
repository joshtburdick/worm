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
cells = cells[ nchar(cells) == 11 ]       # only include terminal cells

cell.col = number.lineage.columns(get.lineage.by.name(lin, "ABplpa"))

m1 = m.unnormalized[ c("ceh-27", "hlh-16", "unc-130", "irx-1"), cells ]

# Plots a lineage in a particular color.
plot.lineage.small = function(gene, main, hue) {
  par(mar=c(0,0.2,1,0.2) + 0.1)

  a = hsv(hue, 1, 0.8 * m.unnormalized[gene,])

  names(a) = lin.node.names

  plot(1,1, xlim=range(cell.col), ylim=c(440,50),
    xaxt="n", yaxt="n", bty="n")
  par(new=FALSE)
  plot.segments.per.cell(a, main, root=lineage.to.plot, times=c(50,440),
      add=TRUE, lwd=5, yaxt="n")
 # axis(2, font=1)
#  mtext("Time (minutes)", side=2, line=2, cex=0.7)

}

pdf("git/sort_paper/plot/figure/figure1.pdf",
  width=8, height=8)

layout(matrix(c(1:5,6,6,6,6,7), nrow=5),
  widths=c(2,1), heights=c(1,1,1,1,1.3))

# plot the trees
plot.lineage.small("ceh-27", "", 0)
plot.lineage.small("hlh-16", "", 1/6)
plot.lineage.small("unc-130", "", 2/6)
plot.lineage.small("irx-1", "", 3/6)

# made-up expression
x = matrix(0, nrow=ncol(m1), ncol=5)
# x[1:20, 1:2] = 1
# x[5:8, 3:4] = 1
# x[34:49, 5] = 0.8
x[ grepl("ABplpaaaa", cells), 1:2 ] = 1
x[ grepl("ABplpaaa", cells), 3 ] = 1
x[ grepl("ABplpapaa", cells), 4:5 ] = 1

x = 0.8*x + rnorm(prod(dim(x)), mean=0.05, sd=0.05)
x[x < 0] = 0
x[x > 1] = 1

par(mar=c(0,0.2,0,0.2) + 0.1)
plot(1,1, type="n", xlim=range(cell.col), ylim=c(0, 5),
  bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
image(x = cell.col[cells], y = c(0:4)+0.5, x, add=TRUE,
  col=hsv(0,0,128:0/128), useRaster=FALSE, bty="n", xaxt="n", yaxt="n")

# color slices of heatmap
color.slice = function(y1, y2, hue) {
  rect(min(cell.col)-0.5, y1, max(cell.col)+0.5, y2,
    col=hsv(hue, 1, 1, 0.1), border=NA)
}
color.slice(0, 2, 0.1)
color.slice(2, 3, 0.3)
color.slice(3, 5, 0.7)

par(mar=c(0,1,0,1))
plot(1,1, type="n", xlim=c(0,4), ylim=c(0.33,4), bty="n",
  main="", xlab="", ylab="", xaxt="n", yaxt="n")

for(i in 0:3) {
  j = i + 0.5
  lines(c(-0.5,j,j), c(j,j,0.2), xpd=NA)
}

x.p = m1 %*% x
x.n = (1-m1) %*% x
x.f = x.p / (x.p + x.n)

x.f = x.f[4:1,]    # to simplify drawing tubes

par(mar=c(0,0.2,0,0.2)+0.1)
plot(1,1, type="n", xlim=c(0,4), ylim=c(0,5), bty="n",
  main="", xlab="", ylab="", xaxt="n", yaxt="n")
image(x=c(0:3)+0.5, y=c(0:4)+0.5, x.f, col=hsv(0,0,128:0/128),
  add=TRUE, useRaster=FALSE, bty="n", xaxt="n", yaxt="n")

# again, color slices
color.slice = function(y1, y2, hue) {
  rect(0, y1, 4, y2,
    col=hsv(hue, 1, 1, 0.1), border=NA)
}
color.slice(0, 2, 0.1)
color.slice(2, 3, 0.3)
color.slice(3, 5, 0.7)

dev.off()

