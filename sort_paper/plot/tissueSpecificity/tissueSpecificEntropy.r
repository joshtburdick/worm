# Measuring how tissue-specific expression is.

source("git/unmix/ept/gamma.r")
source("git/sort_paper/plot/tissueSpecificity/geneGroups.r")

# the read ratios
rpm = as.matrix(read.tsv("git/cluster/readsPerMillion.tsv"))

# just using FACS data, and genes expressed at least somewhat
sf = grep("(HS |RNAi)", colnames(rpm), invert=TRUE, value=TRUE)
rpm = rpm[ , sf ]

rpm = rpm[ apply(log2(rpm+1),1,max) > 0, ]

# measures of how expressed a given gene is
# ??? not obvious which of these to use
max.expr = apply(log2(rpm+1), 1, max)
# mean.expr = apply(log2(rpm+1), 1, mean)

gene.groups = lapply(gene.groups,
  function(x) intersect(x, rownames(rpm)))

# convert to fraction-of-total
# ??? what should the pseudocount be?
# r = (rpm + 3) / apply( (rpm + 3), 1, sum)
r = (rpm + 1e-6) / apply((rpm + 1e-6), 1, sum)

# definition of entropy
entropy = function(x) {
  x1 = ifelse(x == 0, 0, x * log2(x))
  - sum( x1 )
}

tissue.spec.entropy = apply(r, 1, entropy)

pdf("git/sort_paper/plot/tissueSpecificity/tissueSpecificEntropy.pdf",
  width=7.5, height=10)
par(mfrow=c(3,2))

for(group in names(gene.groups)) {
  cat(group)

  # first, stats for all the genes in this group
  x1 = max.expr[ gene.groups[[group]] ]
  h = hist(x1,
    freq=TRUE, main=paste("Expression of", group), cex.main=1.5,
    xlab="Maximum log2(rpm+1)",
    col="lightgrey", breaks=40)

  # moment-matched gamma (we'll use this later)
  g = gamma.mv2s( rbind(m = mean(x1), v = var(x1)) )
  par(new=TRUE)
  plot(0:1000 / 10,
    length(x1) * dgamma(0:1000/10, shape=g["a",], scale=g["b",]),
    xlim=range(x1),
    main="", xlab="", ylab="", xaxt="n", yaxt="n",
    type="l", lwd=2, col="#0000ff80")

  # the tissue-specific entropy for this
  e = tissue.spec.entropy[ gene.groups[[group]] ]
  hist(e, freq=TRUE, main=paste("Entropy of", group), cex.main=1.5,
    xlab="Entropy (bits)", xlim=c(0,log2(ncol(r))),
    col="lightgrey", breaks=30)
  mtext(paste("mean =", signif(mean(e), 2),
    "  median =", signif(median(e), 2),
    "  n =", length(e)))

  # pick a matched set of protein-coding genes
  p = dgamma(max.expr[ gene.groups[["protein-coding"]] ],
      shape=g["a",], scale=g["b",])
  mg = sample(gene.groups[["protein-coding"]],
    length(gene.groups[[group]]), replace=TRUE,
    prob = p)

  # same things, for matched protein-coding genes
  # first, stats for all the genes in this group
  x1 = max.expr[ mg ]
  hist(x1,
    freq=TRUE, main="Expression of matched\nprotein-coding genes",
    cex.main=1.5, xlab="Maximum log2(rpm+1)",
    col="lightgrey", breaks=40)

  # the tissue-specific entropy for this
  e = tissue.spec.entropy[ mg ]
  hist(e, freq=TRUE, main="Entropy of matched\nprotein-coding genes",
    cex.main=1.5, xlab="Entropy (bits)", xlim=c(0,log2(ncol(r))),
    col="lightgrey", breaks=30)
  mtext(paste("mean =", signif(mean(e), 2),
    "  median =", signif(median(e), 2),
    "  n =", length(e)))

  plot.new()
  plot.new()
}

dev.off()

