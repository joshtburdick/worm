# Unmixing, estinating all cells at once.

library(limSolve)
# library(NMF)    not sure this is working

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
m.orig = {

  # add in negative fractions
  neg = 1 - sortMatrix
  rownames(neg) = paste(rownames(sortMatrix), "m", sep=".")
  m = rbind(sortMatrix, neg)

  # add in "all" fraction
  m = rbind(all=1, m)

  m
}

# mapping between rows of the sort matrix, and sequenced fractions
seqAndFractions = read.table("git/unmix/prediction/seqAndFractions.tsv",
  sep="\t", header=TRUE, as.is=TRUE)

# XXX for now, omitting double-sorted genes, as we don't
# have negatives for them (and they may be good for validation)
seqAndFractions = seqAndFractions[ grep("hlh16", seqAndFractions$fraction.name, invert=TRUE) , ]

# subset sort and expression matrices to those genes
m.orig = m.orig[ seqAndFractions$fraction.name , ]
rownames(m.orig) = seqAndFractions$seq.experiment
rna.seq = rna.seq[ , seqAndFractions$seq.experiment ]

# normalize
m = m.orig / apply(m.orig, 1, sum)



# Unmixes a set of genes using an unmixing function.
unmix.all = function(m, x.f, unmix.f, output.file,
    write.interval = 100) {
  x = matrix(nrow=nrow(x.f), ncol=ncol(m))
  rownames(x) = rownames(x.f)
  colnames(x) = colnames(m)

  unmix.result = list(x = x, i = NA)

  for(i in 1:nrow(x.f)) {
    try({
      r = unmix.f(m, as.vector(x.f[i,]))
      unmix.result$x[i,] = r$x
    })

    if (i %% write.interval == 0) {
      unmix.result$i = i
      cat("**** i =", i, "\n")
      save(unmix.result, file=output.file)
    }
  }

  save(unmix.result, file=output.file)
}

# Unmixes one gene using nnls.
# Doesn't include any prior (simply a penalty on the sum.)
unmix.nnls = function(m, x.fraction) {
  w = 1e-3    # weight for penalty for sum of this
  num.cells = dim(m)[2]

  r = nnls(A = rbind(m, w * diag(num.cells)),
    B = c(x.fraction, rep(0, num.cells)), verbose=TRUE)
  x = as.vector(r$X)
  names(x) = colnames(m)
  list(x = x)
}

# Unmixes (all genes) using fcnnls.
# XXX seems to not be running very fast at all.
# Args:
#   m - sort matrix
#   x.f - expression in fractions
#   w - weight to add penalty for length of solution
# Returns: matrix of unmixed expression.
unmix.fcnnls = function(m, x.f, w=1e-2) {

  x = rbind(m, a=1/ncol(m))
  y = rbind(t(x.f), w)
  r = fcnnls(x, y, pseudo=TRUE)

  r
}

# limit to genes with expression at least 100
# reads per million in some fraction
# rna.seq = rna.seq[ apply(rna.seq, 1, max) >= 50 , ]



# unmix.all(m, rna.seq, unmix.nnls,
#   "git/unmix/prediction/20130423/unmix.nnls.Rdata")





