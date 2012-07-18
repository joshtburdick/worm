# Finds cell volume using an embryo imaged with plasma membrane.

image.path = "data/image/tifs/"
output.path = "git/unmix/image/mrf/"

system(paste("mkdir -p ", output.path))

library("rtiff")

r = readTiff(paste(image.path, "red/20080221_OD70-t055-p19.tif", sep="/"))@red
g = readTiff(paste(image.path, "green/20080221_OD70-t055-p19.tif", sep="/"))@green

result.tif = newPixmapRGB(r, g, 0*r)

# the belief states: this is the probability that a pixel is on the membrane
m = array(dim=c(3, nrow(r), ncol(r)))
dimnames(m)[[1]] = c("n", "a", "b")
m["a",,] = 0
m["b",,] = 1
m["n",,] = 1

# Updates the probability that a pixel is on the membrane.
mrf.update = function(m, r, g) {
  nr = nrow(r)
  nc = ncol(r)

  # normalize
  m["a",,] = m["a",,] / m["n",,]
  m["b",,] = m["b",,] / m["n",,]
  m["n",,] = 1

  # create new, and send messages in all four directions
  m1 = 0 * m
  m1[,1:(nr-1),] = m1[,1:(nr-1),] + m[,2:nr,]
  m1[,2:nr,] = m1[,2:nr,] + m[,1:(nr-1),]
  m1[,,1:(nc-1)] = m1[,,1:(nc-1)] + m[,,2:nc]
  m1[,,2:nc] = m1[,,2:nc] + m[,,1:(nc-1)]

  # add in messages from the actual image
  m1["a",,] = m1["a",,] + r
  m1["n",,] = m1["n",,] + r
  m1["b",,] = m1["b",,] + g
  m1["n",,] = m1["b",,] + g

  # normalize
  m1["a",,] = m1["a",,] / m1["n",,]
  m1["b",,] = m1["b",,] / m1["n",,]
  m1["n",,] = 1
  m1
}

run.it = function(m) {
  for(iter in 1:200) {
    cat(iter, "")
    m = mrf.update(m, r, g)
cat(max(m), "")
    result.tif@blue = m["a",,] / m["n",,]

    writeTiff(result.tif, paste(output.path, "mrf.step.", iter, ".tiff", sep=""))
  }

  m
}

# m = run.it(m)

