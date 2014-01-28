# Optimizes the pseudocount, and enrichment cutoff,
# for considering genes "enriched" (or "depleted".)

library("idr")

# note: this re-generates the "readRatios.tsv" table
source("git/cluster/readRatios.r")

cnd1.samples = c("cnd-1 12/14", "cnd-1 1/4", "cnd-1 8/19")
pha4.samples = c("pha-4 12/9", "pha-4 5/9", "pha-4 9/1")

# Computes enrichment with a particular pseudocount added.
compute.enrichment = function(r, samples, pseudocount) {
  r.pos = r[ , paste(samples, "(+)") ]
  r.neg = r[ , paste(samples, "(-)") ]

  enrich = log2( pseudocount + r.pos ) - log2( pseudocount + r.neg )  
  enrich
}

cnd1.enrich = compute.enrichment(readsPerMillion, cnd1.samples, 1)
pha4.enrich = compute.enrichment(readsPerMillion, pha4.samples, 1)

# Gets the correspondence curve.
get.correspondence.curve.1 = function(x1, x2, n=1000) {
  rank.1 = rank(x1)
  rank.2 = rank(x2)
  uv = get.correspondence(rank(x1), rank(x2), (1:n) / (n+1))

  uv
}

# Plots correspondence curves for different numbers of pseudocounts.
# Doesn't set up the coordinate axes.
plot.correspondence = function(r, samples, depletion, pseudocount, col) {

  en = compute.enrichment(r, samples, pseudocount)
  if (depletion) {
    en = -en
  }
  en = en[ apply(en, 1, mean) > 0 , ]

  uv12 = get.correspondence.curve.1(en[,1], en[,2])
  uv13 = get.correspondence.curve.1(en[,1], en[,3])
  uv23 = get.correspondence.curve.1(en[,2], en[,3])

#  xlim = range(c(uv12$dpsi$t, uv13$dpsi$t, uv23$dpsi$t))
#  ylim = range(c(uv12$dpsi$value, uv13$dpsi$value, uv23$dpsi$value))

  lines(uv12$dpsi$smoothed.line, lwd=2, col=col)
  lines(uv13$dpsi$smoothed.line, lwd=2, col=col)
  lines(uv23$dpsi$smoothed.line, lwd=2, col=col)
}

# Plots correspondence curves for a variety of pseudocounts.
plot.correspondence.diff.pseudocounts = function(r, samples, depletion, main) {
  pseudocounts = c(1:5,10)
  colors = hsv(1:6/7, 1, 0.5, alpha=0.5)
  plot(0,0, xlim=c(0,1), ylim=c(0,2), type="n",
    main=main, xlab="fraction of significant peaks", ylab="slope")
  for(i in 1:length(pseudocounts)) {
    plot.correspondence(r, samples, depletion, pseudocounts[i], colors[i])
  }
  legend("topleft", legend=c("Pseudocount", pseudocounts),
    col=c("#00000000", colors), pch=20, cex=0.7)
}

pdf("git/sort_paper/enrichment/idr.pdf", width=10, height=7.5)
par(mfrow=c(2,2))
plot.correspondence.diff.pseudocounts(readsPerMillion,
  pha4.samples, FALSE, "pha-4 enrichment")
plot.correspondence.diff.pseudocounts(readsPerMillion,
  pha4.samples, TRUE, "pha-4 depletion")
plot.correspondence.diff.pseudocounts(readsPerMillion,
  cnd1.samples, FALSE, "cnd-1 enrichment")
plot.correspondence.diff.pseudocounts(readsPerMillion,
  cnd1.samples, TRUE, "cnd-1 depletion")
dev.off()



