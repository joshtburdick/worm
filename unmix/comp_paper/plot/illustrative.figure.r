# Plots some actual expression patterns, and predictions
# in a sub-lineage.

# source("R/unmix/comp_paper/plot/load.unmixed.results.2.r")
source("R/lineage/embryodb.r")
source("git/plot/plot_expr.r")
source("git/plot/smoothing.r")

expr.prediction = unmix.r[["pseudo"]][["expr.cell"]][[30]][,lin.node.names]

# Extends expression so that it's non-decreasing.
# This also extends expression to child nodes.
extend.expression = function(cell, x) {

  x = smooth.dataset(cell, x, cummax)
  max.x.by.cell = c(by(x, cell, max))
  max.x.by.cell["P0"] = 0
  x = pmax(x, max.x.by.cell[parent.of[cell] ], na.rm=TRUE)

  x[is.na(x)] = 0
  max.x.by.cell = c(by(x, cell, max))
  max.x.by.cell["P0"] = 0
  x = pmax(x, max.x.by.cell[parent.of[cell] ], na.rm=TRUE)

  x
}

# Plots part of the C lineage.
# Args:
#   expr - the expression pattern
plot.expr.1 = function(expr, main, max.intensity, max.color) {
  par(mar=c(2,4,1.8,2)+0.1)
  expr[expr > max.intensity] = max.intensity

  c1 = as.vector(col2rgb(max.color) / 255)
  r = rgb(scale.to.unit(expr) * c1[1],
    scale.to.unit(expr) * c1[2],
    scale.to.unit(expr) * c1[3])
  names(r) = names(expr)
  plot.segments.per.cell(r, main, root="C", lwd=3, yaxt="n")
  title(main)
}

# Plots an image.
plot.image = function(series, main, max.intensity, hue) {
  scd = read.embryodb.dat.file(series)

#  scd$expr = scd$local
  expr = smooth.dataset(scd$cell, scd$local, median.smooth(5))
  expr = extend.expression(scd$cell, expr)
  expr[ expr > max.intensity ] = max.intensity
  expr = as.vector(scale.to.unit(expr))
print(range(expr))

  r = data.frame(cell=scd$cell, time.1=scd$time, time.2=scd$time+1,
    col=hsv(hue,1,expr),
#    col  par(mar=c(2,4,1.8,2)+0.1)=rgb(expr,0,0),
    stringsAsFactors=FALSE)
print(r[1:5,])
print(class(r$col))

  par(mar=c(2,4,1.8,2)+0.1)
  plot.segments(r, main, root="C", time=c(0,350), yaxt="n")
}

fig.1 = function() {
  pdf("git/unmix/comp_paper/plot/1. illustration of method.pdf",
    width=8, height=8)
  par(mfcol=c(4,2))
  par(mar=c(2,4,1.8,2)+0.1)

  plot.image("20081016_alr-1_10A2_3_L2", "alr-1", 1200, 0)
  plot.image("20100914_RW10594_nhr-171_L5", "nhr-171", 900, 0.1)
  plot.image("20111003_tlp-1_L3", "tlp-1", 800, 1/3)
  plot.image("20110906_RW10349_lin-39_L1", "lin-39", 500, 0.69)

  plot.image("20110615_pes-1_L4", "pes-1", 400, 0.75)
  plot.new()
  plot.new()
  plot.expr.1(expr.prediction["pes-1",], "pes-1 prediction", 1300,
    hsv(0.75,1,1))

  dev.off()
}

fig.1()


