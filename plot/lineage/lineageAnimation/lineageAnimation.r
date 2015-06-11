# Attempt at drawing a lineage gradually
# (for purposes of displaying with a movie)

# libraries that need to be installed
library(animation)
library(png)
library(tiff)

#### begin configuration

# directory of images, and name, of the movie to draw
movie.dir = "/gpfs/fs0/l/murr/azach/images/20140311_mom-2_L3/"
movie.name = "20140311_mom-2_L3"

# image file of lineage to draw
lineage.image.filename = "20140311_mom-2_L3.png"

# time interval of stacks to include
time.range = c(100, 115)

# corresponding y-coordinates in the lineage image, in pixels
# (where 1 is at top of the image, and 1 is the bottom). This
# will be revealed, linearly, at the same time that the movie
# plays.
y.range.pixels = c(80, 162)

# speed of the movie
frame.rate.in.seconds = 0.5

### end configuration

# Gets a max. projection z-stack.
get.stack = function(dir, name, t) {

  # XXX sizes of images shouldn't be hard-coded
  img = array(0, dim=c(512, 712, 3))

  # loop through the image planes
  for(plane in c(1:67)) {
    red = readTIFF(sprintf("%s/tifR/%s-t%03d-p%02d.tif",
      dir, name, t, plane))
    green = readTIFF(sprintf("%s/tif/%s-t%03d-p%02d.tif",
      dir, name, t, plane))

    img[,,1] = pmax(img[,,1], red)
    img[,,2] = pmax(img[,,2], green)
  }

  # tweak colors
  img[,,1] = 3 * img[,,1]       # red
  img[,,2] = 0.3 * img[,,2]     # green
  img[ img > 1 ] = 1

  as.raster(img)
  #  as.raster(aperm(img, c(2,1,3)))    # if the movie needs rotating
}

lineage.image = readPNG(lineage.image.filename, native=TRUE)

# Draws one image.
draw.images = function() {

  # convert from pixels to "proportion of image"
  y.range = y.range.pixels / nrow(lineage.image)

  for(i in time.range[1]:time.range[2]) {
    cat(i, "")
    img = get.stack(movie.dir, movie.name, i)
    # this gives the relative sizes of the two panels
    layout(matrix(1:2, nrow=2), heights = c(0.8, 0.2))

    # the movie
    par(mar=c(0,0,0,0))
    plot(0,0, xlim=c(0,1), ylim=c(0,1), type="n", bty="n",
         xaxt="n", yaxt="n")
    rasterImage(img, 0.1,0,0.9,1)

    # the lineage
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

  ani.options(interval = frame.rate.in.seconds, nmax = 1000)

  # XXX this isn't working, presumably because ffmpeg is missing
#   saveVideo(draw.images(), ffmpeg="/home/jburdick/gcb/git/plot/lineage/ffmpeg_like",  # XXX
#              ani.width=400, ani.height=300, clean=FALSE)

  saveGIF(draw.images(),
    ani.width=400, ani.height=300)
}

write.movie.with.lineage()

