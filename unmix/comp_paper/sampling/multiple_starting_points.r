# Attempt to pick multiple starting points to
# diagnose convergence (or lack thereof.)

colors1 = hsv(0:9/10, 1, 0.5, alpha=0.95)

# Picks one "overdispersed" starting point, by finding a solution,
# using a prior with a random center (and unit variance.)
# Args:
#   m - the cell-sorting matrix
#   x.f - the expression in each group of cells
#   alpha - amount to "pull" solution towards the edge point
# Returns: a point x, such that m %*% x = x.f, and x >= 0.
overdispersed.start.1 = function(m, x.f, alpha = 0.9) {
  n = ncol(m)

  # first, the "pseudoinverse" starting point ...
  x0 = lsei(A = diag(n), B = rep(0, n),
    E = m, F = x.f, G = diag(n), H = rep(0, n), tol=1e-6)$X

  # ... then, a starting point way out on an edge
  x1 = lsei(A = diag(n), B = runif(n) * max(x.f),
    E = m, F = x.f, G = diag(n), H = rep(0, n), tol=1e-6)$X

  (1-alpha) * x0 + alpha * x1
}

# Reads in one set of "overdispersed" results, and
# summarizes it.
# Args:
#   f - the file to read
# Returns: list of arrays, each of which is one overdispersed result.
read.overdispersed = function(g) {
  unmix.result = NULL
  load(paste("/home/jburdick/gcb/git/unmix/comp_paper/sampling/multiple_restart/", g, ".sampling.Rdata",
    sep="", collapse=""))
  if (!is.null(unmix.result))
    unmix.result$r[[1]]$x.summary
  else
    return(NULL)

}

# Plots just a summary of these.
plot.summary = function(r, gene, ylim) {
  xlim=c(1,1341)

  plot(1,1, xlim=xlim, ylim=ylim, type="n",
    main=paste(gene, ", all iterations"),
    xlab="", ylab="", xaxt="n", yaxt="n")
  for(j in 1:10) {
    y = r[[j]][,"mean",]
    y.mean = apply(y, 2, mean)
    y.sd = apply(y, 2, sd)
    y.lo = y.mean - y.sd
    y.hi = y.mean + y.sd
    par(new=TRUE)
    segments(1:1341, y.lo, 1:1341, y.hi, col=colors1[j], xlim=xlim, ylim=ylim,
      main="", xlab="", ylab="", xaxt="n", yaxt="n")
  }


}

# Plots results from overdispersed starting points.
# Args:
#   r - the statistics about overdispersed starting points
#   g - the name of the gene
# Side effects: writes out the file
plot.overdispersed = function(r, g) {

  ymax = 0    # FIXME
  for(j in 1:10) {
    m = 1.5 * max(r[[j]][,"mean",])
    if (m > ymax)
      ymax <- m
  }

  png(paste(g, "_multiple_starts.png", sep="", collapse=""), width=1200, height=800)
  par(mfrow=c(7,3))
  par(mar=c(2,1,1,0.1)+0.1)
  xlim=c(1,1341)
  ylim=c(0,ymax)

  # plot groups at various times
  for(i in 1:20) {
    plot(1,1, xlim=xlim, ylim=ylim, type="n",
      main=paste(g, ", ", i, "'th group of 5000 iterations"),
      xlab="", ylab="", xaxt="n", yaxt="n")
    for(j in 1:10) {
      y = r[[j]][i,,]
      y.sd = sqrt(y["var",])
      y.lo = y["mean",] - y.sd
      y.hi = y["mean",] + y.sd
      par(new=TRUE)
      segments(1:1341, y.lo, 1:1341, y.hi, col=colors1[j], xlim=xlim, ylim=ylim,
        main="", xlab="", ylab="", xaxt="n", yaxt="n")
    }
  }

  plot.summary(r, g, ylim)

  dev.off()
}

if (FALSE) {
for(g in c("alr-1", "B0310.2", "B0336.3", "C05D10.1", "C08B11.3")) {
  cat(g, "")
  r = read.overdispersed(g)
  plot.overdispersed(r, g)
}
}


