

source("git/unmix/ept/beta.r")
source("git/unmix/ept/gamma.r")



# Moment-matches a beta with a gamma.
beta.gamma.mm = function(a, b) {

  g = gamma.mv2s(beta.s2mv(rbind(a=a,b=b)))

  x = c(0:100) / 100
  ylim = c(0,4)
  plot(x, dbeta(x, a, b), type="l", xlab="", ylab="", ylim=ylim,
    main=paste("Beta(", a, ", ", b, ")", sep=""))
  par(new=TRUE)
  plot(x, dgamma(x, shape=g["a",1], rate=g["b",1]),
    type="l", col="#808080", xlab="", ylab="", ylim=ylim)

}





pdf("git/tex/unmix/gamma/betaGamma.pdf", width=9, height=3)
par(mfrow=c(1,3))
beta.gamma.mm(1,3)
beta.gamma.mm(4,3)
beta.gamma.mm(3,1)
dev.off()

