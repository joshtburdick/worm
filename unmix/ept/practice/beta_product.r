# More fiddling with beta distributions.
# Conjecture: the product of beta distributions should
# look like a beta distribution.

# ... and in fact it does seem to. The scales of the beta
# distributions don't seem to follow a particular pattern,
# however.

source("git/utils.r")
source("git/unmix/ept/beta.r")

# Converts from mean and variance to moments.
mv2moment = function(mv) {
  rbind(x1 = mv["m",],
    x2 = mv["m",] ^ 2 + mv["v",])
}

# ... and back.
moment2mv = function(m) {
  rbind(m = m["x1",],
    v = m["x2",] - m["x1",] ^ 2)
}

# Moment-matching estimate of parameters of the product.
beta.product.moment.match = function(a1, b1, a2, b2) {
  m1 = mv2moment(beta.s2mv(rbind(a=a1, b=b1)))
  m2 = mv2moment(beta.s2mv(rbind(a=a2, b=b2)))
  m = m1 * m2
  beta.mv2s(moment2mv(m))
}

# Plots the product of beta-distributed variables.
plot.beta.product = function(a1, b1, a2, b2, n=1000000) {
  write.status(paste(a1, b1, a2, b2))

  x = rbeta(n, a1, b1)
  y = rbeta(n, a2, b2)
  x1 = c(0:100) / 100

  hist(x, xlim=c(0,1), col="grey", main=paste("(", a1, ",", b1, ")"))
  hist(y, xlim=c(0,1), col="grey", main=paste("(", a2, ",", b2, ")"))

  p = x*y

  # compute parameters of moment-matched estimate
  mm = beta.mv2s(rbind(m = mean(p), v = var(p)))

  hist(p, xlim=c(0,1), col="grey",
    main=paste0("(", round(mm["a",], 3), ", ", round(mm["b",], 3), ")"))

  # plot moment-matched beta which fits the sampled distribution
  par(new=TRUE)
  plot(x1, dbeta(x1, mm["a",], mm["b",]), type="l", lwd=2, col="#00008080",
    xlab="", ylab="", xaxt="n", yaxt="n")
  
  # ...and closed-form approximation
  m = beta.product.moment.match(a1, b1, a2, b2)
  par(new=TRUE)
  plot(x1, dbeta(x1, m["a",], m["b",]), type="l", lwd=2, col="#80000080",
    xlab="", ylab="", xaxt="n", yaxt="n")
  mtext(paste0("predicted = (", round(m["a",], 3), ", ", round(m["b",], 3), ")"),
    cex=0.8, col="#800000b0", line=0)
}


pdf("git/unmix/ept/practice/beta_product.pdf", width=7.5, height=10)
par(mfrow=c(4,3))

# plot.beta.product(1,1, 1,1)
# plot.beta.product(3,1, 1,3)
# plot.beta.product(1,3, 4,4)
# plot.beta.product(2,5, 3,3)

for(x in 1:4)
  for(y in x:4)
    for(z in 1:4)
      for(w in z:4)
        plot.beta.product(x,y,z,w)
dev.off()

foo = beta.product.moment.match(1,2,3,4)
