# Plots of some properties of Dirichlet, conditional on a sum
# of some of them being beta-distributed.

load("git/unmix/ept/practice/dirichlet.beta.samples.Rdata")

a = dirichlet.beta.samples[,c(1:5)]
b = dirichlet.beta.samples[,c(6,7)]
p = dirichlet.beta.samples[,c(8:12)]

s.a = apply(a, 1, sum)
s.b = apply(b, 1, sum)
s.p = apply(p, 1, sum)


a.m = a / s.a
p.m = p / s.p

png("git/unmix/ept/practice/dirichlet_margin_plots.png", width=1200, height=600)
par(mfrow=c(1,2))

# more pointed observation
i = s.b > 14

# the concentration
lim=c(12, 48)
plot(s.a, s.p, pch=20, cex=1.5, col="#00000060", xlim=lim, ylim=lim,
  xlab = "prior concentration", ylab="posterior concentration")
par(new=TRUE)
plot(s.a[i], s.p[i], pch=20, cex=1.5, col="#ff0000b0", xlim=lim, ylim=lim, xlab="", ylab="")
abline(0,1)

# the mean
plot(as.vector(a.m), as.vector(p.m), pch=20, cex=1.5, col="#00000060",xlim=c(0,0.5), ylim=c(0,0.5), xlab="prior mean", ylab="posterior mean")
par(new=TRUE)
plot(as.vector(a.m[i,]), as.vector(p.m[i,]), pch=20, cex=1.5, col="#ff0000b0", xlim=c(0,0.5), ylim=c(0,0.5), xlab="", ylab="")
abline(0,1)

dev.off()

