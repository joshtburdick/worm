# Does unmixing using sampling.

source("unmix_test.r")

source("sampling/cdaCpp.r")

# Does unmixing of one gene.
# Args:
#   num.samples, thinning
# Returns: function which, given a sort matrix and
#   expression in each fraction, returns unmixed expression.
unmix.sampling = function(num.samples, thinning) {

  function(m, x.fraction) {
    num.cells = dim(m)[2]

    # first, find cells estimated to be zero
    x0 = lsei(E=m, F=x.fraction, G=diag(num.cells), H=rep(0, num.cells),
      tol=1e-6, type=2)$X
cat("called lsei\n")
    cl = which(abs(x0) >= 1e-30)
cat(length(cl), "were non-zero\n")

    # do sampling
st = system.time(X1 <- sample.cda(m[,cl], x.fraction, x0[cl], num.samples, thinning))

    # add zero cells back in
    # XXX this wastes memory
    X = matrix(0, nrow=num.samples, ncol=num.cells)
    X[,cl] = X1
    colnames(X) = colnames(m)

    # just storing summary stats...
    x = apply(X[c((num.samples/2):num.samples),], 2, mean)
    s = apply(X[c((num.samples/2):num.samples),], 2, sd)
    names(x) = colnames(m)
    names(s) = colnames(m)
    list(x = x, sd = s, last.sample = X[num.samples,], system.time = st)
  }
}

# Does sampling for each gene.
unmix.all = function(genes) {
  out.dir = "sampling/expr.30/"
  system(paste("mkdir -p ", out.dir))

  for(g in genes) {
    unmix.result = run.unmix.1(expr.cell[g,,drop=FALSE], m.cell, unmix.sampling(20000, 1000), reporters$picked, 30) 
    save(unmix.result, file=paste(out.dir, "/", g, ".Rdata", sep="", collapse=""))
  }
}

# unmix.all(rownames(expr.cell))
# unmix.all(rev(rownames(expr.cell)))

