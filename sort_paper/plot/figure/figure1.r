# Figure 1: showing lineages used for sorting, and the
# corresponding scatterplots.


source("git/unmix/seq/plot/scatterplot.r")

source("git/sort_paper/unmix/sortMatrix.r")
source("git/plot/plot_expr.r")
source("git/plot/label_panel.r")

# font size for panel labels
gp = gpar(fontsize=12, col="black")

# Draws a scatterplot, coloring genes which are far away from
# a threshold.
draw.scatter.diag.threshold.1 =
    function(x, y, x.name, y.name, cutoff = 2) {
  par(mar=c(4,4,3,3)+0.1)
  lim = c(min(x,y), max(x, y))
  main = expression("expression, log"[2](3+coverage))

  colors = ifelse(y-x > cutoff, "#80000080",
    ifelse(y-x < -cutoff, "#00008080", "#00000080"))

  par(font.axis=5)
  plot(x, y, xlab=x.name, ylab=y.name,
    main = main, xlim=lim, ylim=lim, col=colors,
    pch=183, font=5, xaxt="n", yaxt="n",
    cex=0.8, cex.axis=1, cex.lab=1.2, cex.main=1.2)
  axis(1, font=1)
  axis(2, font=1)
  abline(0,1, lwd=2, col="#00000040")
  abline(2,1, lwd=2, col="#80000040")
  abline(-2,1, lwd=2, col="#00008040")
}

pdf("git/sort_paper/plot/figure/figure1.pdf",
  width=10, height=10)

layout(matrix(c(1:6), nrow=3, byrow=TRUE),
  widths=c(2,1), heights=c(1,1,1))

par(mar=c(4,4,3,3)+0.1)
a = rgb(m.unnormalized["pha-4",], 0, 0)
names(a) = lin.node.names
plot.segments.per.cell(a, "pha-4 expression", root="P0", times=c(0,300),
    lwd=1.5, yaxt="n")
axis(2, font=1)
mtext("Time (minutes)", side=2, line=2, cex=0.7)
label.panel("a)", gp=gp)

draw.scatter.diag.threshold.1(r[,"pha-4 9/1 (-)"], r[,"pha-4 9/1 (+)"],
  "pha-4 9/1 (-)", "pha-4 9/1 (+)", 2)
label.panel("d)", gp=gp)

par(mar=c(4,4,3,3)+0.1)
a = rgb(m.unnormalized["ttx-3",], 0, 0)
names(a) = lin.node.names
plot.segments.per.cell(a, "ttx-3 expression", root="P0", times=c(0,300),
    lwd=1.5, yaxt="n")
axis(2, font=1)
mtext("Time (minutes)", side=2, line=2, cex=0.7)
label.panel("b)", gp=gp)

draw.scatter.diag.threshold.1(r[,"ttx-3 (-)"], r[,"ttx-3 (+)"],
  "ttx-3 (-)", "ttx-3 (+)", 2)
label.panel("e)", gp=gp)

par(mar=c(4,4,3,3)+0.1)
a = rgb(m.unnormalized["hlh-16",], m.unnormalized["ceh-6",], 0)
names(a) = lin.node.names
plot.segments.per.cell(a, "hlh-16 and ceh-6 expression", root="ABp", times=c(0,300),
    lwd=2.5, yaxt="n")
axis(2, font=1)
mtext("Time (minutes)", side=2, line=2, cex=0.7)
label.panel("c)", gp=gp)

draw.scatter.diag.threshold.1(
  r[,"ceh-6 (-) hlh-16 (-)"], r[,"ceh-6 (+) hlh-16 (+)"],
  "ceh-6 (-) hlh-16 (-)", "ceh-6 (+) hlh-16 (+)", 2)
label.panel("f)", gp=gp)

dev.off()

