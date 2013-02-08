# Models embryo cell locations using a multivariate normal.

library("mnormt")
library("rgl")

source("R/lineage/embryodb.r")

# Gets the coordinates from a file at some time.
# Args:
#   scd - the scd file
#   t - the timepoint
# Returns: the distance matrix, scaled to microns.
get.points = function(scd, t) {
  a = scd[ scd$time == t , ]
  r = data.frame(x = a$x, y = a$y, z = a$z)
  rownames(r) = a$cell
  r
}

# Scales a matrix, then removes attributes.
scale1 = function(x) {
  r = scale(x)
  attributes(r) <- NULL
  dim(r) = dim(x)
  dimnames(r) = dimnames(x)
  r
}

# Given coordinates, creates a covariance matrix based
# on their distances.
# Args:
#   x - the coordinates
# Returns: a covariance matrix. This probably won't
# usually be s.p.d.
get.dist.covariance = function(x, length.scale=1) {
  d = dist(scale(x))
  s = exp(-( d / length.scale) ^2 )
  kronecker(as.matrix(d)^2, diag(3))
}

# Samples a random embryo.
# Args:
#   si2 - the covariance matrix from 
#   noise - amount of noise variance to add
# Returns: coordinates for points drawn
sample.random.embryo = function(si2, noise) {

  S = si2 + diag(noise, nrow(si2))

  # sample a random vector
  x = as.vector(rmnorm(varcov=S))

  # reshape it
  r = as.data.frame(matrix(x, ncol=3, byrow=TRUE))
  colnames(r) = c("x", "y", "z")
  r
}

# Renders an embryo using rgl.
# Args:
#   a - data.frame or matrix with columns "x", "y", and "z".
#   offset - amount to add to x, y, and z coordinates
# Side effects: draws an embryo.
render.embryo = function(a, offset) {
  spheres3d(a[,"x"] + offset[1], a[,"y"] + offset[2], a[,"z"] + offset[3],
    radius=0.3, color = rainbow(nrow(a)))
}

# Draws an embryo, and some samples from the same distribution.
draw.embryos = function() {
  r = get.points(e1, 140)
  render.embryo(scale(r), c(0, 0, 10))

  si2 = get.dist.covariance(scale1(r))
  for(i in 1:5) {
    r.sampled = scale1(sample.random.embryo(si2, 500))
    render.embryo(r.sampled, c(10 * (i-1), 0, 0))
  }
}

# Plots a histogram of the log-probability for the
# actual assignment, compared to random shufflings
# of the points.
lp.shuffle.hist = function(x, noise, n) {
  x = scale1(x)

  as.column = function(x) as.vector(t(as.matrix(x)))

  si2 = get.dist.covariance(x, 0.1)
  si2 = si2 + diag(noise, nrow(si2))

  l1 = dmnorm(as.column(x), varcov = si2, log = TRUE)
  r = rep(NA, n)
  for(i in 1:n) {
    r[i] = dmnorm(as.column(x[sample(nrow(x)),]),
      varcov = si2, log = TRUE)
  }

  xlim = range(c(r, l1))
  hist(r, xlim = xlim, breaks=50)
  abline(v=l1, col="red", lwd=3)
  print(l1)
}

# Same as above, but uses correlation of distances.
dist.shuffle.hist = function(x, n) {
  x = scale1(x)
  
  d.actual = as.vector(dist(x))

  r = rep(NA, n)
  for(i in 1:n)
    r[i] = cor(d.actual, as.vector(dist(x[sample(nrow(x)),])))

  hist(r, xlim = c(min(r), 1), breaks=100, col="grey")
}

# read a random embryo file
e1 = read.embryodb.dat.file("20111208_JIM108_L1")

r = get.points(e1, 140)
si2 = get.dist.covariance(scale(r))

open3d()
draw.embryos()

if (FALSE) {
  pdf("git/unmix/image/mvnorm_embryo_dist_shuffle.pdf")
  dist.shuffle.hist(r,100000)
  dev.off()
}
