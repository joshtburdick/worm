# Plots expression patterns, per-cell predictions, and
# total predictions within a lineage.

source("R/plot/plot_expr.r")
source("git/plot/label_panel.r")

output.dir = "git/unmix/comp_paper/plot/EP.expr.in.groups/"

expr.cell = as.matrix(read.table("~/gcb/git/unmix/unmix_comp/data/exprCell.tsv",
  sep="\t", header=TRUE, row.names=1, as.is=TRUE))

load("~/gcb/git/unmix/unmix_comp/src/EP.2/lineage.totals.Rdata")

# num.cells = apply(cell.lineage.matrix, 1, sum)

# utility to clip something to be positive
clip.pos = function(x) ifelse(x < 0, 0, x)

system(paste("mkdir -p ", output.dir))

# Plots the uncertainty for a given gene's prediction.
plot.it = function(g) {
  cat(g, "")
  pdf(file=paste(output.dir, "/", g, ".pdf", sep=""), width=11, height=10)
  par(mfrow=c(3,1))
  par( mar=c(2,4,2,1) + 0.1 )
  z = 0 * expr.cell[1,]

  # first, plot the actual expression
  col = rgb(scale.to.unit(expr.cell[g,]), z, z)
  names(col) = names(expr.cell[g,])
  plot.segments.per.cell(col, cell.time.on.off, lwd=1.1)
  label.panel("a)")
  title(paste(g, "expression"))

  x = t(lineage.totals[g,,])

  # then, the per-cell EP prediction mean
  col = rgb(scale.to.unit(x["per.cell.mean",]),
    scale.to.unit(sqrt(clip.pos(x["per.cell.var",]))),
    z)
  names(col) = names(expr.cell[g,])
  plot.segments.per.cell(col, cell.time.on.off, lwd=1.1)
  label.panel("b)")
  title(paste(g, "per-cell mean and SD"))

  # then, lineage mean, colored by amount of uncertainty
  col = rgb(scale.to.unit(x["lineage.mean",]),
    scale.to.unit(sqrt(clip.pos(x["lineage.var",]))),
    z)
  names(col) = names(expr.cell[g,])
  plot.segments.per.cell(col, cell.time.on.off, lwd=1.1)
  label.panel("c)")
  title(paste(g, "per-lineage mean and SD"))

  dev.off()
}

for(g in rownames(expr.cell)) {
  plot.it(g)
}


