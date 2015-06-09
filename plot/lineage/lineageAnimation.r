# Attempt at drawing a lineage gradually
# (for purposes of displaying with a movie)

library(animation)
# library(pixel)
library(grid)
library(png)
library(tiff)

# the movie to draw
movie.dir = "/gpfs/fs0/l/murr/azach/images/20140311_mom-2_L3/"
movie.name = "20140311_mom-2_L3"

# image file of lineage to draw
lineage.image.filename = "image/movie/20140311_mom-2_L3.png"

# time interval to include
time.range = c(50, 200)

# corresponding x-coordinates in the lineage image
# (where 0 is the top of the image, and 1 is the bottom)
y.range = c(0.5, 1)

# Gets a max. projection z-stack.
get.stack = function(dir, name, t) {

  img = array(0, dim=c(512, 712, 3))

  for(plane in c(30:40)) {
    red = readTIFF(sprintf("%s/tifR/%s-t%03d-p%02d.tif",
      dir, name, t, plane))
    green = readTIFF(sprintf("%s/tif/%s-t%03d-p%02d.tif",
      dir, name, t, plane))

    img[,,1] = pmax(img[,,1], red)
    img[,,2] = pmax(img[,,2], green)
  }

  # tweak colors
  img[,,1] = 3 * img[,,1]
  img[,,2] = 0.3 * img[,,2]
  img[ img > 1 ] = 1

  as.raster(img)
#  as.raster(aperm(img, c(2,1,3)))
}

lineage.image = readPNG(lineage.image.filename, native=TRUE)



draw.images = function() {
  for(i in time.range[1]:time.range[2]) {
    img = get.stack(movie.dir, movie.name, i)
    layout(matrix(1:2, nrow=2), heights = c(0.8, 0.2))

    # the movie
    par(mar=c(0,0,0,0))
    plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n", bty="n",
      xaxt="n", yaxt="n")
    rasterImage(img, 0.1,0,0.9,1)

    # the lineage (partially covered)
    par(mar=c(0,0,0,0))
    plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n", bty="n",
      xaxt="n", yaxt="n")
    rasterImage(lineage.image, 0,0,1,1)
    # cover up part of the lineage
    i1 = (i - time.range[1]) / (time.range[2] - time.range[1])
    y1 = i1 * (y.range[2] - y.range[1]) + y.range[1]
    par(new=TRUE)
    rect(0,0,1,1-y1, col="white", lwd=0)
  }
}

# Writes a movie of images, while gradually
# drawing a lineage tree.
write.movie.with.lineage = function() {

  ani.options(interval = 0.1, nmax = 300)
  saveGIF(draw.images(),
    ani.width=1024, ani.height=768)

  # XXX saveGIF() seems to have problems with pathnames;
  # hacking around that
  system(paste("mv animation.gif ",
    "git/plot/lineage/lineageAnimation.gif"))
}


# foo = get.stack(movie.dir, movie.name, 240)
# rasterImage(foo,0,0,1,1)

write.movie.with.lineage()


