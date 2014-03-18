# Optimizes the pseudocount, and enrichment cutoff,
# for considering genes "enriched" (or "depleted".)

# note: this re-generates the "readRatios.tsv" table
source("git/cluster/readRatios.r")

# Computes enrichment, using some pseudocount.
# Args:
#   r - the "reads per million" counts
#   samples - names of samples (without "(+)" or "(-)")
#   pseudocount - the pseudocount of reads per million to add
# Returns: matrix of how much each gene is enriched in each
#   sample (negative numbers indicate depletion)
enrichment = function(r, samples, pseudocount) {
  r.pos = r[ , paste(samples, "(+)") ]
  r.neg = r[ , paste(samples, "(-)") ]

  enrich = log2( pseudocount + r.pos ) - log2( pseudocount + r.neg )  
}

cnd1.samples = c("cnd-1 12/14", "cnd-1 1/4", "cnd-1 8/19")
pha4.samples = c("pha-4 12/9", "pha-4 9/1", "pha-4 5/9")

# Computes how often an enrichment is reproduced.
# Args:
#   x - the enrichment of each gene in one sample
#   y - enrichment of each gene in another sample
# Returns: after ordering genes in y according to x,
#   the cumulative proportion of genes in y which are
#   enriched (in terms of at least having a positive sign)
enrich.reproducibility = function(x, y) {
  stopifnot( length(x) == length(y) )
  i = order(x, decreasing=TRUE)
  y1 = y[ i ]
  s = cumsum( y1 > 0 ) / (1:length(x))
  list(x = x[ i ], y = cumsum( y1 > 0 ) / (1:length(x)) )
}

# some "toy" data sets, which should have the same
# "enrichment reproducibility".
x1 = cbind(c(5:-4), c(5,2,-1,-2,4,5,1,1,1,-1))
x2 = x1[sample(1:10),]

# Computes enrichment reproducibility for all pairs of
# samples.
# Args:
#   en - enrichment of each gene (with one sample per column)
# Returns: list with two matrices (one column per pair of samples)
#   x - matrix of first samples' enrichment
#   y - fraction of things reproducible per pair of samples.
enrich.reproducibility.all.pairs = function(en) {
  n = ncol(en)
  x = NULL
  y = NULL
  for(i in 1:n)
    for(j in 1:n)
      if (i != j) {
        er = enrich.reproducibility(en[,i], en[,j])
        x = cbind(x, er$x)
        y = cbind(y, er$y)
      }
  list(x = x, y = y)
}

plot.enrichment.reproducibility = function(r, samples, depletion, main) {

  en = enrichment(r, samples, 3)
  if (depletion) {
    en = -en
  }
  er = enrich.reproducibility.all.pairs(en)

  plot(1,1, xlim=c(0,5), ylim=c(0.5,1), type="n", main=main,
    xlab="enrichment in first replicate",
    ylab="fraction enriched in second replicate")

  for(j in 1:ncol(er$x)) {
    lines(er$x[,j], er$y[,j], col="#0000ff80")
  }

  abline(v = 2, col="#00000080")
}

pdf("git/sort_paper/enrichment/cutoffOptimize.pdf",
  width=7.5, height=7.5)
par(mfrow=c(2,2))
plot.enrichment.reproducibility(readsPerMillion, pha4.samples,
  FALSE, "pha-4 enrichment")
plot.enrichment.reproducibility(readsPerMillion, pha4.samples,
  TRUE, "pha-4 depletion")
plot.enrichment.reproducibility(readsPerMillion, cnd1.samples,
  FALSE, "cnd-1 enrichment")
plot.enrichment.reproducibility(readsPerMillion, cnd1.samples,
  TRUE, "cnd-1 depletion")
dev.off()


