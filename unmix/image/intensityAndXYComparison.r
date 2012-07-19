# Comparison of image intensities for particular cells, and change in XY position.

source("R/lineage/embryodb.r")

# For a given movie, plots intensity versus time for cells which are
# "on" at least part of the time, and 
plot.intensity.and.loc.diff = function(series.name, intensity.cutoff) {
  a = read.embryodb.dat.file(series.name)

  # compute location differences
  dx = c(NA, diff(a$x))
  dy = c(NA, diff(a$y))
  a$tot.loc.diff = sqrt(dx^2 + dy^2)

  cell.list = unique(a[a$blot>intensity.cutoff,"cell"])

  for(i in sort(cell.list)) {
    plot(a[a$cell==i,"blot"], type="l", xlab="time", ylab="blot intensity",
      main=paste(series.name, i))
    par(new=TRUE)
    plot(a[a$cell==i,"tot.loc.diff"], type="l", col="red",
      xlab="", ylab="", xaxt="n", yaxt="n", ylim=c(0,50))
    axis(side=4, col="red", col.axis="red")
    mtext("movement (pixels)", side=4, col="red", cex=0.75, line=2.5)
  }
}



pdf("git/unmix/image/intensityAndXYComparison.pdf", width=7.5, height=10)
par(mfrow=c(4,2))
par(mar=c(5,4,4,4)+0.1)
plot.intensity.and.loc.diff("20110209_UP2051_mls-2_L2", 7000)
plot.intensity.and.loc.diff("20110726_RW10713_L1", 5000)
plot.intensity.and.loc.diff("20110731_RW10434_L1", 4000)
dev.off()

