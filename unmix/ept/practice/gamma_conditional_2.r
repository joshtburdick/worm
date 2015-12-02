# Attempt at plotting this density function.
# XXX deprecated -- this doesn't seem correct.

source("git/unmix/ept/gamma.r")
source("git/unmix/ept/beta.r")

x0 = c(0:500) / 100


# starting with examples which actually are Dirichlet on the marginal
# two distributions, as gamma parameters
# g1 = rbind(a=c(4,1), b=c(12,3))
# g1 = rbind(a=c(2,3), b=c(10,10))
g1 = rbind(a=c(2,3), b=c(11,13))
g1.n = gamma.s2n(g1)

# ... moment-matched as beta
b1 = beta.mv2s(gamma.s2mv(g1))







# print(gamma.conditional.numerical.1(g1))

b1.flip = rbind(a = b1["b",], b = b1["a",])


bp = (beta.s2n(b1)[,1,drop=FALSE] + beta.s2n(b1.flip)[,2,drop=FALSE]) / 2
bp.s = beta.n2s(bp)

# plot the product
plot(x0, dgamma(x0, g1["a",1], g1["b",1]) * dgamma(1-x0, g1["a",2], g1["b",2]),
  ylab="", type="l")

# ... and what I think the product should be, using beta approximations
par(new=TRUE)
# plot(x0, dbeta(x0, bp.s["a",1], bp.s["b",1]), ylab="", type="l", col="#0000f080", lwd=5)
plot(x0, dbeta(x0, b1["a",1], b1["b",1]) * dbeta(x0, b1.flip["a",2], b1.flip["b",2]),
  ylab="", type="l", col="#0000f080", lwd=5)



# ... and using a gamma approximation
if (TRUE) {
  g1.flip = gamma.mv2s(beta.s2mv(b1.flip))
  par(new=TRUE)
  plot(x0, dgamma(x0, g1["a",1], g1["b",1]) * dgamma(x0, g1.flip["a",2], g1.flip["b",2]),
    ylab="", type="l", col="#f0000080", lwd=5)
}


# print(beta.n2mv(bp))

# a moment-matching approximation (possibly not used)
if (FALSE) {
  mv2mm = function(x) rbind(x1 = x["m",], x2 = x["m",]^2 + x["v",])
  mm2mv = function(x) rbind(m = x["x1",], v = x["x2",] - x["x1",]^2)
  foo = (mv2mm(beta.s2mv(b1[,1,drop=FALSE])) + mv2mm(beta.s2mv(b1.flip[,2,drop=FALSE]))) / 2
  print(mm2mv(foo))
}



# Given some gamma distributions, plots several approximations to the
# distribution of the first (conditional on them summing to 1.)
plot.conditional.approx = function(a) {

  y = gamma.conditional.density(a)(x0)
#  ylim = c(0,0.5*max(y))

  plot(x0, y / max(y),
    xlim=c(0,5), ylim=c(0,1), type="l",
    xlab="", ylab="", lwd=3, col="#00000080")

  # scaling approximation (deprecated)
  if (FALSE) {
  a1 = array(gamma.s2n(a), dim=c(nrow(a), ncol(a), 1))
  dimnames(a1)[[1]] = c("e1","e2")
  a.scaling = gamma.n2s(gamma.conditional.1(a1))
  a.scaling.y = dgamma(x0, a.scaling["a",1], a.scaling["b",1])
  par(new=TRUE)
  plot(x0, a.scaling.y / max(a.scaling.y),
    xlim=c(0,5), ylim=c(0,1), type="l", lwd=3, col="#00800080",
    xlab="", ylab="")
  }

  # beta approximation
  a.beta = gamma.n2s(gamma.conditional.approx.1(gamma.s2n(a)))
  a.beta.y = dgamma(x0, a.beta["a",1], a.beta["b",1])
  par(new=TRUE)
  plot(x0, a.beta.y / max(a.beta.y),
    xlim=c(0,5), ylim=c(0,1), type="l", lwd=3, col="#80000080",
    xlab="", ylab="")
}

# Plots several approximations.
plot.approximations = function() {
  pdf("git/unmix/ept/practice/gamma_conditional_2.pdf",
    title="Approximating gamma, conditional on sum",
    width=7.5, height=10)
  par(mfrow=c(3,2))

  plot.conditional.approx(rbind(a=c(2,3), b=c(10,10)))
  plot.conditional.approx(rbind(a=c(2,3), b=c(30,20)))
  plot.conditional.approx(rbind(a=c(3,2), b=c(20,30)))

  plot.conditional.approx(rbind(a=c(1,10,5), b=c(30,20,16)))
  plot.conditional.approx(rbind(a=c(5,5), b=c(30,15)))
  plot.conditional.approx(rbind(a=c(2,5), b=c(10,20)))
  plot.conditional.approx(rbind(a=c(5,1), b=c(30,30)))

  plot.conditional.approx(rbind(a=c(7,7), b=c(60,40)))
  plot.conditional.approx(rbind(a=c(7,7), b=c(40,55)))

  dev.off()
}

# Plots an approximation of an exponential with a beta.
# Args:
#   b - the scale of the gamma (the shape is 1)
plot.exp.beta.comparison = function(b) {
  x = c(0:1000) / 1000
  y = dgamma(x, shape=1, rate=b)
  ylim = c(0, max(y))
  plot(x, y, ylim=ylim, type="l", lwd=3, col="#50505080")

  p = beta.mv2s(rbind(
    m=trunc.gamma.moment(1, 1, 1, b) / trunc.gamma.moment(0,1,1,b),
    v=trunc.gamma.moment(2, 1, 1, b) / trunc.gamma.moment(0,1,1,b)))
print(p)
  par(new=TRUE)
  plot(x, dbeta(x, p["a",], p["b",]), ylim=ylim, type="l", lwd=3, col="#f0000080")
}

plot.approximations()

