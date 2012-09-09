# Plots a movie as a tree.

source("git/plot/plot_expr.r")
source("R/lineage/embryodb.r")



# Plots a movie.
# Args:
#   series.name - name of the series to plot
#   main - title
#   lwd - with of lines
plot.movie = function(series.name, main, lwd=5) {
  scd = read.embryodb.dat.file(series.name)
  if (is.null(scd)) {
    cat("\nfailed to read series ", series.name, "... skipping\n")
    return
  }

  # only include the cells in the lineage
  scd = scd[scd$cell %in% lin.node.names,]
  r1 = data.frame(cell=scd$cell, time.1=scd$time, time.2 = scd$time+1,
    col=rgb(scale.to.unit(scd$blot), 0, 0), stringsAsFactors=FALSE)

  plot.segments(r1, main, root="P0",
    times=c(0, 350), lwd=lwd)
}

# Plots several movies as PNGs.
# Args:
#   series.names - what series to plot
#   output.dir - where to write output to
#   width, height - size of output images
#   time - times to include
#   lwd - line width
plot.movie.pngs = function(series.names, output.dir,
  width=1000, height=400) {
  system(paste("mkdir -p", output.dir))

  for(s in series.names) {
    cat(s, "")
    png(paste(output.dir, s, ".png", sep=""), width=width, height=height)
    plot.movie(s, s)
    dev.off()
  }
}

plot.movie.pngs(c("20110929_ceh-43_L1"), "~/tmp/movie.pngs/")



