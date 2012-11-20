# Compares error margins in particular sublineages.

cwd = getwd()
setwd("~/gcb/git/unmix/unmix_comp/src/")
source("unmix_test.r")
load("../data/tree_utils.Rdata")
load("EP/lineage.totals.Rdata")
load("EP/lineage.totals.unscaled.Rdata")
setwd(cwd)

num.cells = apply(cell.lineage.matrix, 1, sum)

# utility to clip something to be positive
clip.pos = function(x) ifelse(x < 0, 0, x)

if (FALSE) {

representative.per.gene =
  read.table("data/image/representative_per_gene.tsv",
    header=TRUE, sep="\t", as.is=TRUE)

r = representative.per.gene[
  representative.per.gene$time.series.name %in% rownames(lineage.totals) &
  representative.per.gene$gene %in% rownames(expr.cell) , ]

expr.cell = expr.cell[ r$gene , ]
lineage.totals = lineage.totals[ r$time.series.name , , ]
dimnames(lineage.totals)[[1]] = r$gene

}

# compute actual mean expression in each lineage
actual.mean.expr = t( t(expr.cell %*% t(cell.lineage.matrix)) / num.cells )



# Given a set of lineages, computes statistics about where the
# actual expression falls within those bounds.
plot.expression.bounds = function(lineages, name) {

  zmax = 5

  # prediction mean and standard deviation
if(TRUE) {
  x = 1341 * lineage.totals.unscaled[,lineages,"lineage.mean"]
  y = 1341 * sqrt(clip.pos(lineage.totals.unscaled[,lineages,"lineage.var"]))
  xylim = c(0, max(x, y))
  par(mar=c(5.5,4.5,4,2) + 0.1)
  plot(x, y,
    main=name, xlab="mean of predicted relative expression",
      ylab="SD of predicted relative expression",
    xlim=xylim, ylim=xylim,
    type="p", pch=20, col="#00000080", cex=2,
    cex.axis = 2, cex.lab=2, cex.main=2)
}

  mean.over.sd = lineage.totals[,lineages,"lineage.mean"] /
    sqrt(clip.pos(lineage.totals[,lineages,"lineage.var"]))

#  hist(1 / mean.over.sd, main=paste("sd over mean", name),
#    breaks=100, col="grey", xlim=c(0,40))

  z = (actual.mean.expr[,lineages] - lineage.totals[,lineages,"lineage.mean"]) /
    sqrt(clip.pos(lineage.totals[,lineages,"lineage.var"]))

  z[ z < -zmax ] = -zmax
  z[ z > zmax ] = zmax
  hist(z, main=name, xlab="z", xlim=c(-zmax,zmax),
    breaks=50, col="grey",
    cex.axis = 1.8, cex.lab=1.9, cex.main=1.9)
}

plot.it = function() {
#  pdf("R/unmix/comp_paper/EP/lineage.total.bounds.pdf",
#    width=7, height=10)

#  png("R/unmix/comp_paper/EP/lineage.total.bounds.png", height=1000, width=600)
#  par(mfrow=c(5,2))

  png("git/unmix/comp_paper/plot/lineage.total.bounds.horiz.png", height=900, width=2000)
  par(mfcol=c(2,5))

  plot.expression.bounds(lin.12.cell, "twelve lineages")
  plot.expression.bounds(names(num.cells[num.cells >= 50]),
    "> 50 cells")
  plot.expression.bounds(names(num.cells[num.cells >= 10 & num.cells < 50]),
    "10 - 50 cells")
  plot.expression.bounds(names(num.cells[num.cells >= 5 & num.cells < 10]),
    "5 - 9 cells")
  plot.expression.bounds(names(num.cells[num.cells < 5]),
    "< 5 cells")
  dev.off()
}

plot.it()


