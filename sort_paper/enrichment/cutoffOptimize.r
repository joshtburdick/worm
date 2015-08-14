# Optimizes the pseudocount, and enrichment cutoff,
# for considering genes "enriched" (or "depleted".)

source("git/utils.r")
source("git/sort_paper/plot/experimentRename.r")

readsPerMillion = read.tsv("git/cluster/readsPerMillion.tsv")

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

cnd1.samples = c("cnd-1 8/19", "cnd-1 12/14", "cnd-1 1/4")
pha4.samples = c("pha-4 5/9", "pha-4 9/1", "pha-4 12/9")

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

# Plots enrichment reproducibility, with all the data on one graph.
plot.enrichment.reproducibility.one.graph =
    function(r, samples, depletion, main) {

  en = enrichment(r, samples, 3)
  if (depletion) {
    en = -en
  }
  er = enrich.reproducibility.all.pairs(en)

  par(mar=c(5,5,4,1)+0.1)
  plot(1,1, xlim=c(0,5), ylim=c(0.5,1), type="n", main=main,
    xlab="Enrichment in first replicate",
    ylab="Fraction enriched in\nsecond replicate")
#  mtext("Fraction enriched in\nsecond replicate",
#    side=2, line=2, cex=0.9)

  for(j in 1:ncol(er$x)) {
    lines(er$x[,j], er$y[,j], col="#0000ff80")
  }

  abline(v = 2, col="#00000080")
}

# Plots enrichment reproducibility, broken down by
# individual sample.
plot.enrichment.reproducibility = function(r, samples, depletion, main) {
  en = enrichment(r, samples, 3)
  if (depletion) {
    en = -en
  }

  n = length(samples)
  all.mfe = c()
  for(i in 1:n)
    for(j in 1:n)
      if (i!=j) {

        er = enrich.reproducibility(en[,i], en[,j])

        plot(1,1, xlim=c(0,5), ylim=c(0.5,1), type="n",
          main = if (i==1 && j == 2) main else "",
          cex.main = 1.2,
          xlab=paste(
            ifelse(depletion, "Depletion in", "Enrichment in"),
            prettify.read.ratio.columns(samples[i])),
          ylab=paste(
            ifelse(depletion, "Fraction depleted\nin", "Fraction enriched\nin"), 
            prettify.read.ratio.columns(samples[j])))

        lines(er$x, er$y, col="#0000ffff")
        abline(v = 2, col="#00000080")

        min.fraction.enriched = min( er$y[ er$x >= 2 ] )
        all.mfe = c(all.mfe, min.fraction.enriched)
        mtext(paste("min. fraction",
          ifelse(depletion, "depleted", "enriched"),
          "=", round(min.fraction.enriched, 2)), line=0.3, cex=0.7)
      }

  cat(paste0(main, ": mean min fraction enriched/depleted = ",
    round(mean(all.mfe), 4), "\n"))

}

plot.separate.reproducibility.curves = function() {
pdf("git/sort_paper/enrichment/cutoffOptimize.pdf",
  width=8.5, height=12)
par(mfcol=c(6,4))
par(mar=c(4,5,4,1)+0.1)
plot.enrichment.reproducibility(readsPerMillion, pha4.samples,
  FALSE, "pha-4 enrichment")
plot.enrichment.reproducibility(readsPerMillion, pha4.samples,
  TRUE, "pha-4 depletion")
plot.enrichment.reproducibility(readsPerMillion, cnd1.samples,
  FALSE, "cnd-1 enrichment")
plot.enrichment.reproducibility(readsPerMillion, cnd1.samples,
  TRUE, "cnd-1 depletion")
dev.off()
}

# old version of this, which added a graph to a panel
if (FALSE) {
plot.enrichment.reproducibility.one.graph(readsPerMillion,
  pha4.samples, FALSE, "pha-4 enrichment")
# XXX
# label.panel("c)", gp=gpar(fontsize=12, col="black"))
plot.enrichment.reproducibility.one.graph(readsPerMillion,
  pha4.samples, TRUE, "pha-4 depletion")
plot.enrichment.reproducibility.one.graph(readsPerMillion,
  cnd1.samples, FALSE, "cnd-1 enrichment")
plot.enrichment.reproducibility.one.graph(readsPerMillion,
  cnd1.samples, TRUE, "cnd-1 depletion")
}

plot.separate.reproducibility.curves()

