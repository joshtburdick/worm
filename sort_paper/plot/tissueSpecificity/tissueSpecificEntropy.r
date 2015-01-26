# Measuring how tissue-specific expression is.

source("git/utils.r")
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

# Empirical version of a distribution.
# Args:
#   a - a vector of numbers
#   bucket.size - how to quantize those numbers
# Returns: a function which approximates the distribution of a
empirical.dist = function(a, bucket.size = 0.2) {

  # the probabilities (indexed by bucket)
  p = table(as.character(floor( a / bucket.size )))
  p = p / sum(p)

  function(x) {
    i = as.character(floor( x / bucket.size ))
    ifelse(i %in% names(p), p[i], 0)
  }
}

# Picks a "matched set" of numbers from one set, with distribution
# similar to another set of numbers.
# Args:
#   a - named vector of numbers, whose distribution we want to match
#   b - named vector of numbers, which we will choose from
#   bucket.size - the bucket size
# Returns: vector of n names of elements of b, whose distribution
#   should be similar to that of a
pick.matched = function(a, b, bucket.size = 0.2) {
  f = function(x) as.character(floor( x / bucket.size ))

  buckets = by(names(b), f(b), function(x) c(as.character(x)))

  pick.matched = function(y) {
    sample(buckets[[f(y)]], 1)
  }

  sapply(a, pick.matched)
}

for(group in names(gene.groups)) {
  write.status(group)

  # scatterplot of expression and tissue specificity
  plot(max.expr[ gene.groups[[group]] ],
    tissue.spec.entropy[ gene.groups[[group]] ],
    main = group,
    xlab="Expression (maximum log2(rpm+1))",
    ylab="Entropy (bits)", xaxt="n", yaxt="n",
    pch=183, font=5, col="#0000ff80")
  axis(1)
  axis(2)

  # XXX hack: label some of the "linc"s
  if (group %in% c("ancRNA", "lincRNA", "non-coding")) {
    g = names(max.expr)[ max.expr >= 2 & tissue.spec.entropy <= 4.5 ]
    g = intersect(g, gene.groups[[group]])
    g = setdiff(g, "linc-77")

    text(max.expr[g], tissue.spec.entropy[g],
      labels=paste(" ", g), adj=c(0,0), srt=315,
      cex=ifelse(group=="non-coding", 0.3, 0.6))
  }

  plot.new()

  # first, stats for all the genes in this group
  x1 = max.expr[ gene.groups[[group]] ]
  h = hist(x1,
    freq=TRUE, main=paste("Expression of", group), cex.main=1.4,
    xlab="Maximum log2(rpm+1)",
    col="lightgrey", breaks=40)

  # the tissue-specific entropy for this
  e = tissue.spec.entropy[ gene.groups[[group]] ]
  hist(e, freq=TRUE, main=paste("Entropy of", group), cex.main=1.4,
    xlab="Entropy (bits)", xlim=c(0,log2(ncol(r))),
    col="lightgrey", breaks=30)
  mtext(paste("mean =", signif(mean(e), 2),
    "  median =", signif(median(e), 2),
    "  n =", length(e)))

if (TRUE) {
  # pick a set of protein-coding genes with similar
  # expression levels
  mg = pick.matched(max.expr[ gene.groups[[group]] ],
    max.expr[ gene.groups[["protein-coding"]] ], bucket.size=0.2)

  # same things, for matched protein-coding genes
  # first, stats for all the genes in this group
  x1 = max.expr[ mg ]
  hist(x1,
    freq=TRUE, main="Expression of matched protein-coding genes",
    cex.main=1.4, xlab="Maximum log2(rpm+1)",
    col="lightgrey", breaks=40)

  # the tissue-specific entropy for this
  e = tissue.spec.entropy[ mg ]
  hist(e, freq=TRUE, main="Entropy of matched protein-coding genes",
    cex.main=1.4, xlab="Entropy (bits)", xlim=c(0,log2(ncol(r))),
    col="lightgrey", breaks=30)
  mtext(paste("mean =", signif(mean(e), 2),
    "  median =", signif(median(e), 2),
    "  n =", length(e)))
}
}

dev.off()

# write out table for lincRNAs
lg = gene.groups[["lincRNA"]]
linc.rna.ts = data.frame(max.expr = max.expr[lg],
  tissue.spec.entropy = tissue.spec.entropy[lg])
linc.rna.ts = linc.rna.ts[
  order(-linc.rna.ts$max.expr, -linc.rna.ts$tissue.spec.entropy), ]
write.tsv(linc.rna.ts, "git/sort_paper/plot/tissueSpecificity/lincRNA tissue specificity.tsv")


