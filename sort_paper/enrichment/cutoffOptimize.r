# Optimizes the pseudocount, and enrichment cutoff,
# for considering genes "enriched" (or "depleted".)

# note: this re-generates the "readRatios.tsv" table
# source("git/cluster/readRatios.r")

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
pha4.samples = c("pha-4 12/9", "pha-4 5/9", "pha-4 9/1")

# Computes how often an enrichment is reproduced.
# Args:
#   x - the enrichment of each gene in one sample
#   y - enrichment of each gene in another sample
# Returns: after ordering genes in y according to x,
#   the cumulative proportion of genes in y which are
#   enriched (in terms of at least having a positive sign)
enrich.reproducibility = function(x, y) {
  stopifnot( length(x) == length(y) )
  y1 = y[ order(x, decreasing=TRUE) ]
  s = cumsum( y1 > 0 ) / (1:length(x))
  s
}

# some "toy" data sets, which should have the same
# "enrichment reproducibility".
x1 = cbind(c(5:-4), c(5,2,-1,-2,4,5,1,1,1,-1))
x2 = x1[sample(1:10),]
stopifnot(all( enrich.reproducibility(x1[,1], x1[,2]) == 
  enrich.reproducibility(x2[,1], x2[,2]) ))

# Computes enrichment reproducibility for all pairs of
# samples.
# Args:
#   x - enrichment of each gene (with one sample per column)
# Returns: matrix of reproducibility, with one column
#   per pair of samples.
enrich.reproducibility.all.pairs = function(x) {
  n = ncol(x)
  r = NULL
  for(i in 1:n)
    for(j in 1:n)
      if (i != j) {
        r = cbind(r, enrich.reproducibility(x[,i], x[,j]))
      }
  r
}

plot.enrichment.reproducibility = function(r, samples, depletion, main) {
  # only plot this many genes
  num.genes = 4000

  pseudocount = c(1,5,10,100)
  color = hsv(0:3/4, 1, 1, alpha=0.5)

  plot(1,1, xlim=c(num.genes, 1), ylim=c(0.6,1), type="n", main=main,
    xlab="rank of gene in first replicate", ylab="fraction replicated")

  for(i in 1:length(pseudocount)) {
    en = enrichment(r, samples, pseudocount[i])
    if (depletion) {
      en = -en
    }
    er = enrich.reproducibility.all.pairs(en)
    er = er[ 1:num.genes , ]

    for(j in 1:ncol(er)) {
      lines(1:num.genes, er[,j], col=color[i])
    }

  }

  legend("bottomright", legend = c("Pseudocount", pseudocount),
    col = c("#00000000", color), cex=0.8, pch=20)
}

pdf("git/sort_paper/enrichment/cutoffOptimize.pdf",
  width=7.5, height=10)
par(mfrow=c(4,1))
plot.enrichment.reproducibility(readsPerMillion, pha4.samples,
  FALSE, "pha-4 enrichment")
plot.enrichment.reproducibility(readsPerMillion, pha4.samples,
  TRUE, "pha-4 depletion")
plot.enrichment.reproducibility(readsPerMillion, cnd1.samples,
  FALSE, "cnd-1 enrichment")
plot.enrichment.reproducibility(readsPerMillion, cnd1.samples,
  TRUE, "cnd-1 depletion")
dev.off()


