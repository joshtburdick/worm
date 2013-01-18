# Attempt to pick multiple starting points to
# diagnose convergence (or lack thereof.)

expr.cell = as.matrix(read.table("~/gcb/git/unmix/unmix_comp/data/exprCell.tsv",
  sep="\t", header=TRUE, row.names=1, as.is=TRUE))

colors1 = hsv(0:4/5, 1, 0.5, alpha=0.6)

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
  x0 = lsei(A = diag(n), B = rep(0, n), type=2,
    E = m, F = x.f, G = diag(n), H = rep(0, n), tol=1e-6)$X

  # ... then, a starting point way out on an edge
  x1 = lsei(A = diag(n), B = runif(n) * max(x.f), type=2,
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
  load(paste("git/unmix/comp_paper/sampling/multiple_restart/", g, ".Rdata",
    sep="", collapse=""))
  if (!is.null(unmix.result))
    unmix.result$r[[1]]$r
  else
    return(NULL)
}

# Plots just a summary of these. (Not currently used.)
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
  output.dir = "git/unmix/comp_paper/sampling/multiple_start_plots/"
  system(paste("mkdir -p", output.dir))

  num.starts = length(r)
  num.windows = dim(r[[1]]$x.summary)[1]

cat("num.starts = ", num.starts, "   num.windows = ", num.windows, "\n")

  # first, find maximum y
  ymax = 0
  for(j in 1:num.starts) {
    m = max(r[[j]]$x.summary[,"mean",] + 2 * sqrt(r[[j]]$x.summary[,"var",]))
    if (m > ymax)
      ymax <- m
  }

  pdf(paste(output.dir, "/", g, "_multiple_starts.pdf", sep="", collapse=""), width=7.5, height=10)
  par(mfrow=c(4,1))
  par(mar=c(2,1,1,0.1)+0.1)
  xlim=c(1,1341)
  ylim=c(0,ymax)

  # plot groups at various times
  for(i in 1:num.windows) {
    plot(1,1, xlim=xlim, ylim=ylim, type="n",
      main=paste(g, ", ", i, "'th group of iterations"),
      xlab="", ylab="", xaxt="n", yaxt="n")
    for(j in 1:num.starts) {
      y = r[[j]]$x.summary[i,,]
      y.sd = sqrt(y["var",])
      y.lo = y["mean",] - y.sd
      y.hi = y["mean",] + y.sd
      par(new=TRUE)
      segments(1:1341, pmax(0, y.lo), 1:1341, y.hi, col=colors1[j], xlim=xlim, ylim=ylim,
        main="", xlab="", ylab="", xaxt="n", yaxt="n")
    }
  }

#  plot.summary(r, g, ylim)

  dev.off()
}

# r = read.overdispersed("alr-1")

if (FALSE) {
for(g in dimnames(expr.cell)[[1]]) {
  cat(g, "")
  r = read.overdispersed(g)
  if (!is.null(r)) {
    plot.overdispersed(r, g)
  }
}
}

