# Trying to compute error bounds.

library("corpcor")

source("R/unmix/comp_paper/plot/load.unmixed.results.2.r")


# expression, and predictions, normalized to sum to one for each gene
# ??? rename this?
rel.expr.cell = expr.cell / apply(expr.cell, 1, sum)
rel.x.d = rel.expr.cell %*% t(m1) %*% pseudoinverse(t(m1))  # x.d / apply(x.d, 1, sum)
rel.x.d.err = rel.x.d - rel.expr.cell


# plot(as.vector(rel.x.d), abs(as.vector(rel.x.d.err)))

# FIXME: use a density plot?
plot.it = function(a, b, a.name, b.name) {
  xy.lim = c(min(a,b), max(a,b))
  print(paste(a.name, b.name))
  plot(as.vector(a), as.vector(b), pch=20, cex=0.1,
    xlim=xy.lim, ylim=xy.lim, xlab=a.name, ylab=b.name,
    main = paste("cor. =", round(cor(as.vector(a), as.vector(b)), 3)))
}


# ??? do these make sense?
per.cell.res = t(m1) %*% pseudoinverse(t(m1))

# Given a matrix, constructs a matrix with just the diagonal elements
# along the diagonal. (Only includes the upper-right corner)
# XXX annoyingly slow
matrix.on.diagonal = function(n) {

  # compute coordinate sets each matrix
  a.index = NULL
  b.index = NULL
  for(s in 1:n) {
if (s %% 50 == 0) cat(s, "")
    a.index = rbind(a.index,
      cbind(rep((1:s), 2), rep((n-s)+(1:s), 2))) 
    b.index = rbind(b.index,
      cbind(n+1-rep(rep(s, s), 2), (n-s) + c(2*(1:s)-1, 2*(1:s))))
  }

#  print(list(a.index=a.index, b.index=b.index))
  function(a) {
    b = matrix(nrow=n, ncol=2*n)
    b[b.index] = a[a.index]
    b
  }
}

# XXX shouldn't have hard-coded numbers
# m.rotate = matrix.on.diagonal(1341)

# Plots actual and simulated expression patterns, along with matrices
# which may or may not be indicative of the structure of errors.
plot.actual.predicted.and.resolution =
  function(x.actual, x.predicted, res.matrix, output.path) {

  library("grid")
#  m.rotate = matrix.on.diagonal(dim(res.matrix)[1])

  for(i in 1:dim(x.actual)[1]) {
cat(rownames(x.actual)[i], "")
    png(paste(output.path, "/", rownames(x.actual)[i], ".png", sep=""),
      width=1500, height=1200)
    par(mfrow=c(4,1))

    c1 = rgb(scale.to.unit(x.actual[i,]), 0, 0)
    names(c1) = colnames(x.actual)
    plot.segments.per.cell(c1)

    c1 = rgb(scale.to.unit(x.predicted[i,]), 0, 0)
    names(c1) = colnames(x.actual)
    plot.segments.per.cell(c1)

    plot(0, 0, xlim=c(0,1), ylim=c(0,1), type="n",
      xaxt="n", yaxt="n", xlab="", ylab="")
    x.p = x.predicted[i,]
    b = 1.5 * (x.p %o% x.p) * res.matrix
    b = b / sd(as.vector(b))
    b[b<0] = 0
    b[b>1] = 1
    b = 1-b
    rasterImage(m.rotate(b), 0, 0, 1, 1, interpolate=FALSE)

    plot(0, 0, xlim=c(0,1), ylim=c(0,1), type="n",
      xaxt="n", yaxt="n", xlab="", ylab="")
    a = 0.5 * res.matrix / sd(as.vector(res.matrix))
    a[a<0] = 0
    a[a>1] = 1
    a = 1-a
    rasterImage(m.rotate(a), 0, 0, 1, 1, interpolate=FALSE)

    dev.off()
  }
}

if (FALSE) {
  # XXX this should include row and column names!
  expr.cell.d = 
  rownames(expr.cell.d) = rownames(m)[-1]
  colnames(expr.cell.d) = colnames(m)
  expr.cell.d = expr.cell.d[rownames(expr.cell),]


  plot.actual.predicted.and.resolution(expr.cell, expr.cell.d, per.cell.res, "out/resolution")
}

