# Trying soft k-means clustering.

library(animation)
library(tiff)

# img.file = "20130717_JIM220_L1-t030-p35.tif"
img.file = "small.tif"

img = readTIFF(img.file)

size = c(x = ncol(img), y = nrow(img))


# Initial centers.
init.m = function(n)
  cbind(x = size["x"] * runif(n),
    y = size["y"] * runif(n),
    v = rgamma(n, shape=1, scale=100))



# Computes (log of) a normal density on a grid.
# (Numbering is 1-based.)
# Args: a - vector of parameters for the normal
# Returns: a matrix, with the log of that normal density.
grid.normal.density = function(a) {
  r = matrix(0, nrow=size["y"], ncol=size["x"])
  r = r + dnorm(1:size["y"], mean=a["y"], sd=sqrt(a["v"]), log=TRUE)
  r = t( t(r) + dnorm(1:size["x"], mean=a["x"], sd=sqrt(a["v"]), log=TRUE))
  r
}

# Computes mean and variance, using moment-matching.
# Args:
#   x - the thing to compute moments of
# Returns: vector of mean and variance of row coordinate.
grid.moments = function(x) {
  n = nrow(x)
  x0 = sum(x)
  x1 = sum(x * c(1:n))
  x2 = sum(x * c(1:n)^2)
  m = x1 / x0

  c(mean = m, var = x2 / x0 - m^2)
}

# Updates once.
# Args:
#   x - the data, as a matrix
#   m - the model parameters
# Returns: list, containing
#   r - array of "responsibilities"
#   m - updated model parameters
ep.update = function(x, m) {

  # compute "responsibilities"
  r = array(dim=c(size["y"], size["x"], nrow(m)))
  for(i in 1:nrow(m))
    r[,,i] = grid.normal.density(m[i,])

  # scale to avoid underflow, convert from log, and normalize
  r = r - as.vector(apply(r, c(1,2), max))
  r = exp(r)
  r = r / as.vector(apply(r, c(1,2), sum))

  # update model parameters
  for(i in 1:nrow(m)) {
    a = (r[,,i] * x)
    x1 = grid.moments(t(a))
    y1 = grid.moments(a)
    m[i,] = c(x = x1["mean"], y = y1["mean"],
      v = (x1["var"] + y1["var"]) / 2)
  }

  list(r = r, m = m)
}

# Plots one timestep.
plot.1 = function(img, r, m, iter) {
  n = dim(r)[3]

  h = c(1:n) / n

  par(mar=c(0,0,0,0))   # XXX
  plot(-1, -1, xlim=c(0,1), ylim=c(0,1))
  par(new=TRUE)
  image(z=img, col=hsv(0,0, c(0:255)/255), useRaster=TRUE)

#  for(i in c(1:n)) {
#    par(new=TRUE)    
#    image(z=r[,,i], col=hsv(h[i], 1, c(0:255)/255, 0.1),
#      zlim=c(0,1), useRaster=TRUE)
#  }

  par(new=TRUE)
  s = sqrt(m[,"v"])
  rect((m[,"y"] - s) / size["y"],
    (m[,"x"] - s) / size["x"],
    (m[,"y"] + s) / size["y"],
    (m[,"x"] + s) / size["x"],
    border = hsv(h, 1, 1, 0.5))
#  points(m[,"x"], m[,"y"], col=hsv(h, 1, 0.5, 0.8))
  legend("bottomright", legend=iter, bty="n", text.col="white")
}

# Does several updates, drawing as it goes.
plot.run = function(img, k, num.steps) {
  m = ep.update(img, init.m(k))

  for(iter in c(1:num.steps)) {
# print(m$m)
    cat(iter, "")
    plot.1(img, m$r, m$m, iter)
    m = ep.update(img, m$m)
  }

}

saveGIF(plot.run(img, 15, 100),
  ani.width=size["y"], ani.height=size["x"],
  interval=0.2)



