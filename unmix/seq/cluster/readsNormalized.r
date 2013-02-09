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

embryo.timeseries = as.matrix(read.table(
  "git/unmix/seq/timing/embryo.timeseries.tsv.gz",
  sep="\t", header=TRUE, row.names=1, as.is=TRUE))
embryo.timeseries =
  t( t(embryo.timeseries) / (apply(embryo.timeseries, 2, sum) / 1e6) )

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
    "cnd-1 8/19", "cnd-1 12/14", "cnd-1 1/4",
    "F21D5.9", "mir-57", "mls-2", "pal-1",
    "pha-4 9/1", "pha-4 12/9", "ttx-3", "unc-130")
  for(g in pos.neg.genes) {
    x = get.average(paste(g, "(+)")) - get.average(paste(g, "(-)"))
    r = cbind(r, g = x)
  }
  colnames(r) = pos.neg.genes
  r
}

# Gets the ratios for genes which lack a matched control.
get.pos.singlet.expression = function(genes) {
  r = NULL
  for(g in genes) {
    x = get.average(g) - r.control[,"singlets"]
    r = cbind(r, g = x)
  }
  colnames(r) = genes
  r
}

# Gets heat-shock compared to time-matched controls.
get.hs = function(g) {
  name.1hr = paste("HS", g, "1hr")
  name.6hr = paste("HS", g, "6hr")

  r = cbind(get.average(name.1hr) - get.average("HS N2 1hr"),
    get.average(name.6hr) - get.average("HS N2 6hr"))

  colnames(r) = c(name.1hr, name.6hr)
  r
}


r.pos.neg = get.pos.neg.expression()

if (FALSE) {
r.pos.singlet = get.pos.singlet.expression(
  c("hlh-16 (+)", "irx-1 (+)",
    "ceh-6 (-) hlh-16 (-)",
    "ceh-6 (-) hlh-16 (+)",
    "ceh-6 (+) hlh-16 (-)",
    "ceh-6 (+) hlh-16 (+)"))
}

singlet.average =
  (get.average("cnd-1 singlets") + get.average("pha-4 singlets")) / 2

# Cases in which it's not obvious what to use as the "negative" sort.
r.no.neg = cbind(
  "hlh-16" = get.average("hlh-16 (+)") -
    get.average("ceh-6 (-) hlh-16 (-)"),
  "irx-1" = get.average("irx-1 (+)") - singlet.average,
  "ceh-6 (+) hlh-16 (+)" = get.average("ceh-6 (+) hlh-16 (+)") -
    get.average("ceh-6 (-) hlh-16 (-)"),
  "ceh-6 (+) hlh-16 (-)" = get.average("ceh-6 (+) hlh-16 (-)") -
    get.average("ceh-6 (-) hlh-16 (-)"),
  "ceh-6 (-) hlh-16 (+)" = get.average("ceh-6 (-) hlh-16 (+)") -
    get.average("ceh-6 (-) hlh-16 (-)"),
  "cnd-1 singlets" = get.average("cnd-1 singlets")
    - get.average("cnd-1 ungated"),
  "pha-4 singlets" = get.average("pha-4 singlets")
    - get.average("pha-4 ungated"))

r.rnai = cbind(
  "RNAi lit-1" = get.average("RNAi lit-1") -
    get.average("RNAi ges-1"),
  "RNAi pop-1" = get.average("RNAi pop-1") -
    get.average("RNAi ges-1"))

r.hs = cbind(get.hs("ceh-32"),
  get.hs("ceh-36"),
  get.hs("elt-1"),
  get.hs("pes-1"),
  get.hs("pha-4"))

r.embryo = t( scale( t(log2(0.1 + embryo.timeseries)), center=TRUE, scale=FALSE) )
r.embryo = rbind(r.embryo, "Y74C9A.3"=NA)   # XXX hack

r.normalized = cbind(r.pos.neg, r.no.neg, r.rnai, r.hs, r.embryo[rownames(r.pos.neg),])
r.normalized[ is.na(r.normalized) ] = 0

write.table(round(r.normalized, 3),
  file="git/unmix/seq/cluster/readsNormalized.tsv",
  sep="\t", row.names=TRUE, col.names=NA)

