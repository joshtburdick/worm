# Just the flow-sorting and timeseries reads.

source("git/utils.r")
source("git/data/biomart_utils.r")

options(stringsAsFactors = FALSE)

wb.gene = useMart("WS220", "wormbase_gene")

experimentNames = read.table("git/unmix/seq/quant/experimentNames.tsv",
  sep="\t", header=TRUE, row.names=1)

count.path = "git/unmix/seq/quant/readsPerMillion"

# reads from FACS experiments
r1 = as.matrix(
  read.table(paste(count.path, "readsPerMillion_Murray_050912.tsv", sep="/"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

r2 = as.matrix(
  read.table(paste(count.path, "readsPerMillion_092812.tsv", sep="/"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

r3 = as.matrix(
  read.table(paste(count.path, "readsPerMillion_20110922.tsv", sep="/"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

r = log2(1 + cbind(r1, r2, r3))
r = r[ , order(colnames(r)) ]

colnames(r) = experimentNames[ colnames(r), "name" ]

# Gets the ratios for all genes which have positive and negative samples.
get.pos.neg.expression = function() {
  a = NULL
  pos.neg.genes = c("ceh-26", "ceh-27", "ceh-36", "ceh-6",
    "cnd-1 8/19", "cnd-1 12/14", "cnd-1 1/4",
    "F21D5.9", "mir-57", "mls-2", "pal-1",
    "pha-4 5/9", "pha-4 9/1", "pha-4 12/9", "ttx-3", "unc-130")
  for(g in pos.neg.genes) {
    x = r[ , paste(g, "(+)") ] - r[ , paste(g, "(-)") ]
    a = cbind(a, g = x)
  }
  colnames(a) = pos.neg.genes
  a
}

singlet.average =
  (r[,"cnd-1 singlets"] + r[,"pha-4 singlets"]) / 2

r.facs = get.pos.neg.expression()
r.facs = cbind(r.facs,
  "hlh-16" = r[,"hlh-16 (+)"] - singlet.average,
  "irx-1" = r[,"irx-1 (+)"] - singlet.average,
  "ceh-6 (+) hlh-16 (+)" =
    r[,"ceh-6 (+) hlh-16 (+)"] - r[,"ceh-6 (-) hlh-16 (-)"],
  "ceh-6 (+) hlh-16 (-)" =
    r[,"ceh-6 (+) hlh-16 (-)"] - r[,"ceh-6 (-) hlh-16 (-)"],
  "ceh-6 (-) hlh-16 (+)" =
    r[,"ceh-6 (-) hlh-16 (+)"] - r[,"ceh-6 (-) hlh-16 (-)"],
  "cnd-1 singlets" = r[,"cnd-1 singlets"] - r[,"cnd-1 ungated"],
  "pha-4 singlets" = r[,"pha-4 singlets"] - r[,"pha-4 ungated"])

# timeseries
r.ts = read.tsv("git/unmix/seq/timing/deconvolved_embryo_ts.tsv")
r.ts.scaled = t( scale( t(log2(1 + r.ts)), center=TRUE, scale=FALSE) )

# more row name munging
gene.names.1 =
  mart.convert(wb.gene, "sequence_name", "public_name")(rownames(r.facs))
rownames(r.facs)[ !is.na(gene.names.1) ] =
  gene.names.1[ !is.na(gene.names.1) ]
gene.names.1 =
  mart.convert(wb.gene, "sequence_name", "public_name")(rownames(r.ts.scaled))
rownames(r.ts.scaled)[ !is.na(gene.names.1) ] =
  gene.names.1[ !is.na(gene.names.1) ]

g = intersect(rownames(r.facs), rownames(r.ts.scaled))

r.ft = cbind(r.facs[g,], r.ts.scaled[g,])

write.tsv(round(r.ft,4), "git/unmix/seq/cluster/readsFACSandTS.tsv")

