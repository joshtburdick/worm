# Attempt at modelling the membrane boundary as a function.


image.path = "data/image/tifs/"
output.path = "git/unmix/image/"

system(paste("mkdir -p ", output.path))

library("rtiff")

# Pads a matrix with 0's.
pad.matrix = function(a, w) {
  r = matrix(0, nrow = nrow(a) + 2*w+1, ncol = ncol(a) + 2*w+1)
  r[ ((w+1):(nrow(a)+w)), ((w+1):(ncol(a)+w)) ] = a
  r
}

r = readTiff(paste(image.path, "red/20080221_OD70-t055-p19.tif", sep="/"))@red
g = readTiff(paste(image.path, "green/20080221_OD70-t055-p19.tif", sep="/"))@green

result.tif = newPixmapRGB(r, g, 0*r)

# Convolves an image with a filter.
# XXX note that this is inefficient for Gaussian blurring.
convolve.filter = function(x, w) {
  s = nrow(w)
  s2 = s/2
  x1 = pad.matrix(x, s)

  r = 0 * x
  for(i in 1:nrow(x)) {
    cat(i, "")
    for(j in 1:ncol(x))
      r[i,j] = sum( x1[ ((i+s2):(i+s2+s-1)),((j+s2):(j+s2+s-1)) ] * w )
  }
  r
}

# A Gaussian filter.
gaussian.filter = function(s, sd) {
  w = matrix(nrow=s, ncol=s)
  s2 = s/2 + 0.5
  for(i in 1:s)
    for(j in 1:s)
      w[i,j] = dnorm(sqrt((i-s2)^2 + (j-s2)^2), mean=0, sd=sd)
  w
}

# Filters an image.
filter.image = function(img, output.file) {
  w = gaussian.filter(40, 20)
  membrane = convolve.filter(img, w)
  membrane[membrane < 0] = 0
  membrane = membrane / max(membrane)
  result.tif@blue = membrane
  writeTiff(result.tif, output.file)
}

filter.image(r, paste(output.path, "gaussian_blur_red.tif"))
filter.image(g, paste(output.path, "gaussian_blur_green.tif"))

