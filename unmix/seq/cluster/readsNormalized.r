# Normalizes read counts, in a hopefully sensible way.
# - if there are pos. and neg. sorted fractions, then pos. - neg.
# - if there are only pos. sorted fractions, then positive - singlets
# - similarly for RNAi: RNAi - singlets
# - for HS: compared to time-matched control

options(stringsAsFactors = FALSE)

experimentNames = read.table("git/unmix/seq/quant/experimentNames.tsv",
  sep="\t", header=TRUE, row.names=1)

count.path = "git/unmix/seq/quant/readsPerMillion"

r1 = as.matrix(
  read.table(paste(count.path, "readsPerMillion_Murray_050912.tsv", sep="/"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

r2 = as.matrix(
  read.table(paste(count.path, "readsPerMillion_092812.tsv", sep="/"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

r = log2(0.1 + cbind(r1, r2))      # ??? add something smaller than 1 here?
r = r[ , order(colnames(r)) ]

# Gets the average for everything some name.
# Args: name - name of an experiment
# Returns: average of all samples with that "name" entry
#   in the experimentNames table.
get.average = function(name) {
  s = rownames(experimentNames[experimentNames$name==name,])
  apply(r[,s,drop=FALSE], 1, mean)
}

# Averages for the controls.
r.control = cbind(
  singlets = apply(r[,c("cndm1_singlets","pham4_singlets")], 1, mean),
  ungated = apply(r[,c("cndm1_ungated","pham4_ungated")], 1, mean))

# Gets the ratios for all genes which have positive and negative samples.
get.pos.neg.expression = function() {
  r = NULL
  pos.neg.genes = c("ceh-26", "ceh-27", "ceh-36", "ceh-6",
    "cnd-1", "F21D5.9", "mir-57", "mls-2", "pal-1",
    "pha-4", "ttx-3", "unc-130")
  for(g in pos.neg.genes) {
    x = get.average(paste(g, "(+)")) - get.average(paste(g, "(-)"))
    r = cbind(r, g = x)
  }
  colnames(r) = pos.neg.genes
  r
}

# Gets the ratios for genes which lack a matched control.
get.pos.singlet.expression = function() {
  r = NULL
  genes = c("hlh-16 (+)", "irx-1 (+)",
    "ceh-6 (-) hlh-16 (-)",
    "ceh-6 (-) hlh-16 (+)",
    "ceh-6 (+) hlh-16 (-)",
    "ceh-6 (+) hlh-16 (+)")
  for(g in genes) {
    x = get.average(g) - r.control[,"singlets"]
    r = cbind(r, g = x)
  }
  colnames(r) = genes
  r
}

r.pos.neg = get.pos.neg.expression()
r.pos.singlet = get.pos.singlet.expression()



write.table(round(r.pos.neg, 3),
  file="git/unmix/seq/cluster/readsNormalized.tsv",
  sep="\t", row.names=TRUE, col.names=NA)

