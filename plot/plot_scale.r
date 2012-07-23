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
  xcol="green", ycol="red", xlim, ylim, xlab, ylab,
  cex) {

  a = matrix(rep(0:8, 9) / 8, nrow=9)
  col = matrix(rgb(as.vector(1-a), as.vector(t(a)), 0, maxColorValue=1), nrow=9)
print(col)

  grid.raster(x=x, y=y, width=width, height=height,
    col)

}


