# Unmixing predictions, including positives and negatives.

source("git/unmix/ept/approx_region.r")

# the expression data, as reads per million
rna.seq = {
  rna.seq.20120509 = as.matrix(read.table(
    "git/unmix/seq/quant/readsPerMillion_Murray_050912.tsv",
    sep="\t", header=TRUE, row.names=1, as.is=TRUE))
  rna.seq.20120928 = as.matrix(read.table(
    "git/unmix/seq/quant/readsPerMillion_092812.tsv",
    sep="\t", header=TRUE, row.names=1, as.is=TRUE))
  cbind(rna.seq.20120509, rna.seq.20120928)
}

# the sort matrix, somewhat processed
sortMatrix = as.matrix(read.table("git/unmix/image/sort/sortMatrix.tsv",
  sep="\t", header=TRUE, row.names=1, as.is=TRUE))
m = {

  # add in negative fractions
  neg = 1 - sortMatrix
  rownames(neg) = paste(rownames(sortMatrix), "m", sep=".")
  m = rbind(sortMatrix, neg)

  # tack on some double-sorted fractions
  m = rbind(m,
    ceh6m.hlh16m = m["ceh-6.m",] * m["hlh-16.m",],
    ceh6m.hlh16p = m["ceh-6.m",] * m["hlh-16",],
    ceh6p.hlh16m = m["ceh-6",] * m["hlh-16.m",],
    ceh6p.hlh16p = m["ceh-6",] * m["hlh-16",])

  # add in "all" fraction
  m = rbind(all=1, m)

  # normalize
  m = m / apply(m, 1, sum)

  m
}

# mapping between rows of the sort matrix, and sequenced fractions
seqAndFractions = read.table("git/unmix/prediction/seqAndFractions.tsv",
  sep="\t", header=TRUE, as.is=TRUE)

# subset sort and expression matrices to those genes
m = m[ seqAndFractions$fraction.name , ]
rownames(m) = seqAndFractions$seq.experiment
rna.seq = rna.seq[ , seqAndFractions$seq.experiment ]

# select genes which are expressed at least 100 rpm
# in some fraction (which admittedly is conservative)
max.expr = apply(rna.seq, 1, max)
rna.seq = rna.seq[ max.expr >= 100 , ]


# Does unmixing using EP.
unmix.ep = function(m, x.f) {

  # for now, normalizing this
  scale = max(x.f)
  x.f = x.f / scale

#  hack which improves convergence
  eps = 1e-4
  x.f = as.vector(x.f + m %*% (rep(eps, ncol(m))))   

  ep = approx.region(m, x.f, eps + x.f,
    converge.tolerance = 1e-9, prior.var = 100, max.iters = 100)

  x = ep$m - rep(eps, ncol(m))

  list(x = x, t = ep$t, update.stats = ep$update.stats,
    scale = scale, eps = eps,
    reporters = rownames(m), x.f = x.f)
}



# foo = unmix.ep(m, rna.seq[2,])
# plot(foo$x)

plot.it = function(g) {
  r = unmix.ep(m, rna.seq[g,] / 1e6)
  plot(r$x[lin.node.names], main=g, type="h")
}


