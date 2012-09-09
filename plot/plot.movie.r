# Plots a movie as a tree.

source("git/plot/plot_expr.r")
source("git/plot/smoothing.r")
source("R/lineage/embryodb.r")

# Plots a movie.
# Args:
#   series.name - name of the series to plot
#   main - title
#   lwd - with of lines
plot.movie = function(series.name, main, lwd=5,
  intensity.bounds=c(0,2000)) {
  scd = read.embryodb.dat.file(series.name)
  if (is.null(scd)) {
    cat("\nfailed to read series ", series.name, "... skipping\n")
    return
  }

  # only include the cells in the lineage
  scd = scd[scd$cell %in% lin.node.names,]

  x.smooth = smooth.dataset(scd$cell, scd$blot, median.smooth(9))

  r1 = data.frame(cell=scd$cell, time.1=scd$time, time.2 = scd$time+1,
    col=rgb(scale.interval.to.unit(x.smooth, intensity.bounds), 0, 0),
    stringsAsFactors=FALSE)

  plot.segments(r1, main, root="P0",
    times=c(0, 350), lwd=lwd)
}

# Plots several movies as PNGs.
# Args:
#   movie.list - data frame with columns
#     series - name of the series
#     name - title for the graph
#     intensity.lo, intensity.hi
#   series.names - what series to plot
#   output.dir - where to write output to
#   width, height - size of output images
#   time - times to include
#   lwd - line width
plot.movie.pngs = function(movie.list, output.dir,
  width=1200, height=400) {
  system(paste("mkdir -p", output.dir))

  for(i in 1:nrow(movie.list)) {
    a = movie.list[i,]
    cat(a$series, a$name, "\n")
    png(paste(output.dir, a$series, ".png", sep=""), width=width, height=height)
    plot.movie(a$series, a$name, lwd=4,
      intensity.bounds=c(a$intensity.lo, a$intensity.hi))
    dev.off()
  }
}

movie.list = read.table("git/plot/movie_list.tsv",
  sep="\t", header=TRUE, as.is=TRUE)

plot.movie.pngs(movie.list, "~/tmp/movie.pngs/")



