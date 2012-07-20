# Finds cell volume using an embryo imaged with a plasma membrane 

image.path = "data/image/tifs/"
output.path = "git/unmix/image/mrf/"

system(paste("mkdir -p ", output.path))

library("rtiff")

r = readTiff(paste(image.path, "red/20080221_OD70-t055-p19.tif", sep="/"))@red
g = readTiff(paste(image.path, "green/20080221_OD70-t055-p19.tif", sep="/"))@green

result.tif = newPixmapRGB(r, g, 0*r)

# the belief states: this is the probability that a pixel is on the membrane
x = array(0, dim=c(3, nrow(r), ncol(r)))
dimnames(x)[[1]] = c("n", "a", "b")

# messages from the actual data
m.image = x
m.image["a",,] = r
m.image["b",,] = g
m.image["n",,] = r + g

# The state of the model, consisting of beliefs and messages in all
# four directions. We assume [1,1] is in the upper-left.
m = list(x = m.image, m.up = x, m.down = x, m.left = x, m.right = x)
# ??? initializing things to "off"
m$x["a",,] = 0
m$x["b",,] = 1
m$x["n",,] = 1

# Normalizes a message.
normalize.message = function(m) {
  m["a",,] = m["a",,] / m["n",,]
  m["b",,] = m["b",,] / m["b",,]
  m["n",,] = 1
  m[ is.na(m) ] = 0
  m
}

# Updates the messages and x.
mrf.update = function(m) {
  nr = nrow(r)
  nc = ncol(r)
  
  # update messages
  m$m.up[,1:nr-1,] = normalize.message( m$x[,2:nr,] - m$m.up[,1:nr-1,] )
  m$m.down[,2:nr,] = normalize.message( m$x[,1:nr-1,] - m$m.down[,2:nr,] )
  m$m.left[,,1:nc-1] = normalize.message( m$x[,,2:nc] - m$m.left[,,1:nc-1] )
  m$m.right[,,2:nc] = normalize.message( m$x[,,1:nc-1] - m$m.right[,,2:nc] )

  # update belief
  m$x = m.image + m$m.up + m$m.down + m$m.left + m$m.right
  m
}

run.it = function(m) {
  for(iter in 1:100) {
    cat(iter, "")
    m = mrf.update(m)
    result.tif@blue = m$x["a",,] / m$x["n",,]
    writeTiff(result.tif, paste(output.path, "mrf.step.", iter, ".tiff", sep=""))
  }

  m
}

# m = run.it(m)

