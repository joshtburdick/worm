# Measuring how tissue-specific expression is by
# the total of the absolute value of the enrichment.

source("git/sort_paper/plot/tissueSpecificity/geneGroups.r")

# the read ratios
rpm = as.matrix(read.tsv("git/cluster/readsPerMillion.tsv"))

# just using FACS data; also include a pseudocount
sf = grep("(HS |RNAi)", colnames(rpm), invert=TRUE, value=TRUE)

# convert to fraction-of-total
r = (rpm + 3) / apply( (rpm + 3), 1, sum)
# r = (rpm + 1e-3) / apply( (rpm + 1e-3), 1, sum)

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
  e = tissue.spec.entropy[gene.groups[[group]] ]
  e = e[ !is.na(e) ]
  hist(e, freq=TRUE, main=group, cex.main=1.5,
    xlab="entropy", xlim=c(0,7),
    col="darkgrey", breaks=30)
  mtext(paste("mean =", signif(mean(e), 2),
    "  median =", signif(median(e), 2),
    "  n =", length(e)))
}
dev.off()

