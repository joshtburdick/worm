# An attempt at a generalization of the Dirichlet.

# XXX This does indeed seem to be a strict generalization
# of a Dirichlet. However, it doesn't seem general enough
# to represent something like, e.g., constraining
# three variables x1, x2, and x3 s.t.
#   x1 + x2 + x3 = 1
#   x1 + x2 = 0.5





# Utility converting from three positive numbers
# to coordinates in two dimensions, in a triangle.
ternary.coords = function(x) {
  # normalize rows of x
  x = x / apply(x, 1, sum)

  # based on, um, Wikipedia
  cbind(x = x[,2] + x[,3] / 2,
    y = sqrt(3/4) * x[,3])
}

# Plots a density as a function of a
# proportion of three variables.
# Args: list with components
#   a1, b1, a2, b2, a3, b3 - parameters of the gamma
#     distributions (in terms of shape and scale)
#   x - matrix with three columns
# Side effects: plots proportions of columns of
#   x, as a triangle.
plot.proportion.3 = function(a, col="#00000040") {

  x1 = ternary.coords(a$x)

#  smoothScatter(x1[,"x"], x1[,"y"], bty="c",
#    nrpoints=0,
#    xaxt="n", yaxt="n", xlab="", ylab="",
#    xlim=c(0,1), ylim=c(0,1))

  plot(x1[,"x"], x1[,"y"],
    bty="n", xaxt="n", yaxt="n", xlab="", ylab="",
    xlim=c(-0.1,1.1), ylim=c(-0.1,1.1),
    pch=183, font=5, cex=1, col=col)
  lines(c(0,1,0.5,0), c(0,0,sqrt(3)/2,0), col="#00000050")

  a1 = round(a$a1, 2)
  b1 = round(a$b1, 2)
  a2 = round(a$a2, 2)
  b2 = round(a$b2, 2)
  a3 = round(a$a3, 2)
  b3 = round(a$b3, 2)

  # label distributions
  text(0, -0.02, paste0("Gamma(", a1, ",", b1, ")"), adj=c(0.2,1), cex=0.9)
  text(1, -0.02, paste0("Gamma(", a2, ",", b2, ")"), adj=c(0.8,1), cex=0.9)
  text(1/2, sqrt(3)/2 + 0.03, paste0("Gamma(", a3, ",", b3, ")"), adj=c(0.5,0), cex=0.9)
}



rgamma3 = function(a1, b1, a2, b2, a3, b3, n=1000) {
  list(a1=a1, b1=b1, a2=a2, b2=b2, a3=a3, b3=b3,
    x = cbind(x1 = rgamma(n, shape=a1, scale=b1),
      x2 = rgamma(n, shape=a2, scale=b2),
      x3 = rgamma(n, shape=a3, scale=b3)))
}

pdf("git/unmix/ept/practice/gamma_conditional_5.pdf",
  width=6, height=8)
par(mfrow=c(4,3))
par(mar=c(0,0,0,0)+0.1)

# plain old Dirichlet
plot.proportion.3(rgamma3(1,1, 1,1, 1,1))
plot.proportion.3(rgamma3(1,1, 1,1, 5,1))
plot.proportion.3(rgamma3(1,1, 5,1, 5,1))

plot.proportion.3(rgamma3(4,1, 5,1, 6,1))
plot.proportion.3(rgamma3(10,1, 10,1, 10,1))
plot.proportion.3(rgamma3(0.1,1, 0.1,1, 0.1,1))

# plot.proportion.3(rgamma3(1,1, 1,1, 3,3))
# plot.proportion.3(rgamma3(1,1, 3,3, 3,3))
# plot.proportion.3(rgamma3(1,3, 3,1, 3,1))

# plot.proportion.3(rgamma3(1,3, 1,3, 3,1))
# plot.proportion.3(rgamma3(1,5, 5,1, 5,1))
# plot.proportion.3(rgamma3(1,1, 5,5, 5,5))

# plot.proportion.3(rgamma3(5,1, 5,5, 5,5))
# plot.proportion.3(rgamma3(1,10, 5,5, 5,5))
plot.proportion.3(rgamma3(1,10, 5,5, 5,5))

# plot.proportion.3(rgamma3(1,1, 10,1, 10,1))
plot.proportion.3(rgamma3(1,10, 10,1, 10,1))
plot.proportion.3(rgamma3(1,100, 50,1, 50,1))

# example of multiplying the gamma distributions together
plot.proportion.3(rgamma3(1,100, 100,0.8, 100,0.2),
  col="#c0000020")
plot.proportion.3(rgamma3(100,0.3, 100,0.7, 1,100),
  col="#00c00020")
plot.proportion.3(rgamma3(
  1+100-1, 1/(1/100+1/0.3),
  100+100-1, 1/(1/0.8+1/0.7),
  100+1-1, 1/(1/0.2+1/100)), col="#c0c00020")

# plot.proportion.3(rgamma3(10000000,0.0001, 10,100, 10,100))

dev.off()
