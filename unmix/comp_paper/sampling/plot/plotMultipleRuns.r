
source("git/plot/plot_expr.r")

source("git/unmix/comp_paper/sampling/convergenceStats.r")

expr.cell = as.matrix(read.table(
    "~/gcb/git/unmix/unmix_comp/data/exprCell.tsv",
  sep="\t", header=TRUE, row.names=1, as.is=TRUE))

colors1 = hsv(0:4/5, 1, 0.5, alpha=0.6)

# Reads in one set of "overdispersed" results, and
# summarizes it.
# Args:
#   f - the file to read
# Returns: list of arrays, each of which is one overdispersed result.
read.overdispersed = function(g) {
  filename = paste(
    "git/unmix/comp_paper/sampling/multiple_restart/", g, ".Rdata",
    sep="", collapse="")
  if (file.exists(filename)) {
    load(filename)
    lapply(unmix.result$r[[1]]$r, function (x) x$x.summary)
  }
  else
    NULL
}

# Computes variance of combined variables.
# Args:
#   means, vars - stats about groups of numbers, with
#     one row per group, and one column per estimand
#   window.size - size of the groups of numbers
# Returns: variance of all the numbers.
#   This isn't quite accurate, but seems fairly close;
# I'm assuming the difference is rounding error.
combined.var = function(means, vars, window.size) {
  stopifnot( nrow(means) == nrow(vars) )
  n1 = window.size * nrow(means)

  v1 = vars * ((window.size-1) / window.size)

  # compute second moment
  x2 = v1 + means^2

  as.vector(apply( x2 - means^2 , 2, mean) * (n1 / (n1 - 1)))
}

# Test of combined.var .
combined.var.test = function() {
  x = matrix(runif(100000, max=c(1,3,5,2,4)), ncol=100, byrow=TRUE)
  window.size = 100

  m = matrix(nrow=10, ncol=100)
  v = matrix(nrow=10, ncol=100)

  for(i in 1:10) {
    i1 = window.size * (i-1) + 1
    i2 = window.size * i
    m[i,] = apply(x[c(i1:i2),], 2, mean) 
    v[i,] = apply(x[c(i1:i2),], 2, var)
  }

  list(x = x,
    var = combined.var(m, v, window.size),
    actual.var = apply(x, 2, var))
}

# Given an array of statistics about windows of sampling,
# computes the combined statistics about all of them.
# Args:
#   a - array of sampling statistics
#   window.size - number of samples aggregated in each window
# Returns: array of statistics. Note that the
#   median is not the actual median, but rather
#   a Tukey-style median-of-medians estimate.
aggregate.sampling.stats.1 = function(r, window.size) {
  rbind(
    mean = apply(r[,"mean",], 2, mean),
    var = as.vector(combined.var(r[,"mean",], r[,"var",], window.size)),
    median = apply(r[,"median",], 2, median),
    min = apply(r[,"min",], 2, min),
    max = apply(r[,"max",], 2, max))
}

# Given several runs, computes statistics on a portion of them.
# Args:
#   r - list of arrays of sampling results
#   window.size - number of samples aggregated in each window
#   w - number of windows to include (note that 
#     equally many burn-in iterations will be skipped, and so
#     the array should have at least twice this many windows.
# Returns: array with dimensions:
#   run - which run this was
#   stat - which statistic (see aggregate.sampling.stats.1())
#   estimand - which thing being estimated (in this case, cell)
aggregate.sampling.stats = function(r, window.size, w) {
  n = length(r)

  a = array(dim=c(n, 5, dim(r[[1]])[3]),
    dimnames=list(run=c(1:n),
      stat=dimnames(r[[1]])[[2]],
      cell=dimnames(r[[1]])[[3]]))
print(dim(a))
print(dim(r[[1]]))
  for(i in 1:length(r)) {
    r1 = r[[i]][c((w+1):(2*w)),,]
    a[i,,] = aggregate.sampling.stats.1(r1, window.size)
  }

  a
}

# Plots results from overdispersed starting points.
# Args:
#   output.dir - where to write the output to
#   r - the statistics about overdispersed starting points
#   g - the name of the gene
# Side effects: writes out the file
plot.overdispersed = function(output.dir, r, g) {
  # this is groups of 2 million iterations (following
  # an equal number of burn-in iterations)
  num.iters = c(1,5,10,25)

  system(paste("mkdir -p", output.dir))

  num.starts = length(r)
  num.windows = dim(r[[1]])[1]

# cat("num.starts = ", num.starts, "   num.windows = ", num.windows, "\n")

  # first, find maximum y
  ymax = 0
  for(j in 1:num.starts) {
    m = max(r[[j]][,"mean",] + 1.0 * sqrt(r[[j]][,"var",]))
    if (m > ymax)
      ymax <- m
  }

  pdf(paste(output.dir, "/", g, "_multiple_starts.pdf", sep="", collapse=""), width=7.5, height=10)   # width=1920, height=1200)
  par(mfrow=c(4,1))
  par(mar=c(4,2,1,0.1)+0.1)
  xlim=c(1,1341)
  ylim=c(0,ymax)

  # plot groups at various times
  for(i in num.iters) {
    plot(1,1, xlim=xlim, ylim=ylim, type="n",
      main=paste(g, ", ", 2*i, "million samples"),
      xlab="", ylab="", xaxt="n", yaxt="n")
    int.n.small = setdiff(int.n, c("Ea", "Ep", "Da", "Dp", "MSpa", "MSpp"))
    axis(1, at=cell.to.column[int.n.small],
      labels=int.n.small,
      las=2, cex.axis=0.6)
    mtext("expression", side=2, cex=0.65)

    for(j in 1:num.starts) {
      y1 = r[[j]][c(i:(2*i)),,,drop=FALSE]
      y = aggregate.sampling.stats.1(y1, 1e6)

      y.sd = sqrt(y["var",])
      y.lo = y["mean",] - 2 * y.sd
      y.hi = y["mean",] + 2 * y.sd
#      segments(1:1341, pmax(0, y.lo), 1:1341, y.hi, col=colors1[j], xlim=xlim, ylim=ylim,
#        main="", xlab="", ylab="", xaxt="n", yaxt="n")
      par(new=TRUE)
      plot(1:1341, pmax(0, y.lo), cex=0.4, pch=20, col="#ff000080",
        xlim=xlim, ylim=ylim, main="",
        xlab="", ylab="", xaxt="n", yaxt="n")
      par(new=TRUE)
      plot(1:1341, pmax(0, y.hi), cex=0.4, pch=20, col="#0000ff80",
        xlim=xlim, ylim=ylim, main="",
        xlab="", ylab="", xaxt="n", yaxt="n")
    }

    # plot actual expression
    par(new=TRUE)
    plot(1:1341, pmax(0, expr.cell[g,]), cex=0.4, pch=20, col="#000000a0",
      xlim=xlim, ylim=ylim, main="",
      xlab="", ylab="", xaxt="n", yaxt="n")
  }


  dev.off()
}

r = read.overdispersed("alr-1")

if (FALSE) {
output.dir = "git/unmix/comp_paper/sampling/multiple_start_plots/"
for(g in c("alr-1")) {    # dimnames(expr.cell)[[1]]) {
  cat(g, "")
  r = read.overdispersed(g)
  if (!is.null(r)) {
    plot.overdispersed(output.dir, r, g)
  }
}
}

s = aggregate.sampling.stats(r, 2e6, 25) 
foo = scale.reduction(s[,"mean",], s[,"var",], 2e6)

