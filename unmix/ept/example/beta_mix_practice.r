# Mixing beta distributions.

output.dir = "git/unmix/ept/example/beta_mix_practice/"

source("git/unmix/ept/approx_region_beta.r")

# Plots a histogram of a beta distribution,
# with a moment-matched beta.
plot.beta.hist = function(x, label) {

  # we're using parameters based on the data (which may introduce
  # slight inaccuracies)
  r = beta.moments.to.params(rbind(mu=mean(x), si2=var(x)))
  a = r["a",1]
  b = r["b",1]
  hist(x, breaks=100, col="grey", xlim=c(0,1),
    main=paste(label, ": a =", round(a,4), "   b =", round(b,4)))
  par(new=TRUE)
  plot(1:999/1000, dbeta(1:999/1000, a,b), type="l",
    main="", xlim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n",
    col="#804080", lwd=2)
}

# ??? this doesn't seem right.
# OK, maybe it is...
plot.beta.mix = function(name, a1,b1, a2,b2, am,bm, n=1e7) {
  system(paste("mkdir -p ", output.dir))

  png(paste(output.dir, "/", name, ".png", sep=""), width=800, height=700)
  par(mfrow=c(2,2))

  xa = rbeta(n, a1, b1)
  xb = rbeta(n, a2, b2)
  xm = rbeta(n, am, bm)
  plot.beta.hist(xa, "A")
  plot.beta.hist(xb, "B")
  plot.beta.hist(xm, "M")

  x = xm * xa + (1-xm) * xb
  plot.beta.hist(x, "mix")

  dev.off()
}

plot.beta.mix("1", 3,2, 2.5,1.3, 1,2.5)
plot.beta.mix("2", 3,2, 2.5,1.3, 100,250)
plot.beta.mix("3", 1,2, 2,1, 3,4)
plot.beta.mix("4", 1,2, 2,1, 6,8)
plot.beta.mix("5", 1,2, 2,1, 1,1)
plot.beta.mix("6", 1,2, 2,1, 10,10)



if (FALSE) {
png("beta_mix.png", width=1200, height=600)
par(mfrow=c(2,1))
hist(x1, col="grey", breaks=100,
  main = paste("mean =", mean(x1), "  var =", var(x1)))
hist(x2, col="grey", breaks=100,
  main = paste("mean =", mean(x2), "  var =", var(x2)))
dev.off()
}


