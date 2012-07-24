# Plots a scale on a graph.

library(grid)



# Draws a 2-D scale on a graph. For now, this uses
# green on the x-axis, and red on the y-axis.
# Args:
#   x, y, width, height - specify location of the scale
#   xlim, ylim - limits on the x- and y- axes
#   xlab, ylab - labels
#   cex - text size for labels
plot.red.green.scale = function(x, y, width, height,
  xcol="green", ycol="red", xmax, ymax, xlab, ylab,
  cex) {

  pushViewport(viewport(x=x, y=y, width=width, height=height, name="scale"))

  a = matrix(rep(0:8, 9) / 8, nrow=9)
  col = matrix(rgb(as.vector(1-a), as.vector(t(a)), 0, maxColorValue=1), nrow=9)

  grid.raster(x=0.5, y=0.5, width=1, height=1, col)
  grid.xaxis(c(0,1), label=c(0,xmax))
  grid.yaxis(c(0,1), label=c(0,ymax))
  grid.text(xlab, x=0.5, y=0, vjust="top") 
  grid.text(ylab, x=0, y=0.5, hjust="right") 

  popViewport()
}


