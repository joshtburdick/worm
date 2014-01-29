# Essentially like readsFACSandTS.r, except that it
# masks missing data (cases where there aren't enough reads,
# and the linc- and anr- transcripts for which we don't have
# numbers from the timeseries data.)

source("git/utils.r")
# source("git/data/biomart_utils.r")
# wb.gene = useMart("WS220", "wormbase_gene")
# ens.ce = useMart("ensembl", "celegans_gene_ensembl")

options(stringsAsFactors = FALSE)

experimentNames = read.table("git/unmix/seq/quant/experimentNames.tsv",
  sep="\t", header=TRUE, row.names=1)

count.path = "git/unmix/seq/quant/readsPerMillion/WS220_20140111/"

# pseudocount to add to read counts
pseudocount = 3

# reads from FACS experiments
r1 = as.matrix(
  read.table(paste(count.path, "050912.tsv", sep="/"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

r2 = as.matrix(
  read.table(paste(count.path, "092812.tsv", sep="/"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

r3 = as.matrix(
  read.table(paste(count.path, "20110922.tsv", sep="/"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

# reads per million
readsPerMillion = cbind(r1, r2, r3)
colnames(readsPerMillion) =
  experimentNames[ colnames(readsPerMillion), "name" ]
readsPerMillion = readsPerMillion[ , order(colnames(readsPerMillion)) ]

write.tsv(round(readsPerMillion,3), "git/cluster/readsPerMillion.tsv")

# log-transformed reads per million, with a pseudocount added
r = log2(pseudocount + readsPerMillion)
r = r[ , order(colnames(r)) ]

# Gets the ratios for all genes which have positive and negative samples.
get.pos.neg.expression = function() {
  a = NULL
  pos.neg.genes = c("ceh-26", "ceh-27", "ceh-36", "ceh-6",
    "cnd-1 8/19", "cnd-1 12/14", "cnd-1 1/4",
    "F21D5.9", "mir-57", "mls-2", "pal-1",
    "pha-4 5/9", "pha-4 9/1", "pha-4 12/9", "ttx-3", "unc-130")
  for(g in pos.neg.genes) {
    pos = r[ , paste(g, "(+)") ]
    neg = r[ , paste(g, "(-)") ]

    x = pos - neg

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
r.ts.scaled = t( scale( t(log2(pseudocount + r.ts)), center=TRUE, scale=FALSE) )

# row name matching; possibly deprecated, as WormMart WS220
# seems to be down
if (FALSE) {
  gene.names.1 =
    mart.convert(wb.gene, "sequence_name", "public_name")(rownames(r.facs))
  rownames(r.facs)[ !is.na(gene.names.1) ] =
    gene.names.1[ !is.na(gene.names.1) ]
  gene.names.1 =
    mart.convert(wb.gene, "sequence_name", "public_name")(rownames(r.ts.scaled))
  rownames(r.ts.scaled)[ !is.na(gene.names.1) ] =
    gene.names.1[ !is.na(gene.names.1) ]
}

# alternative way of matching, using WormBase annotation
# XXX using more recent annotation
gene.ids = read.table(gzfile("data/wormbase/c_elegans.PRJNA13758.WS240.xrefs.txt.gz"), sep="\t")[,c(2,3,4)]
colnames(gene.ids) = c("wb.gene", "gene", "transcript")
# omit "unknown gene name" cases
gene.ids = gene.ids[ gene.ids$gene != "." , ]

# renames a vector of gene names
rename.gene.names = function(a) {
  r = rownames(a)
  i = intersect(r, gene.ids[,"transcript"])
  r[ match(i, r) ] =
    gene.ids[ match(i, gene.ids[,"transcript"]), "gene" ]

  a = a[ !duplicated(r) , ]
  r = r[ !duplicated(r) ]
  rownames(a) = r
  a
}

r.facs = rename.gene.names(r.facs)
r.ts = rename.gene.names(r.ts)
r.ts.scaled = rename.gene.names(r.ts.scaled)

g.nc = grep("^(anr-|linc-)", rownames(r.facs), value=TRUE)

g = union(union(rownames(r.facs), rownames(r.ts.scaled)), g.nc)

r.ft = cbind(as.data.frame(r.facs)[g,], as.data.frame(r.ts.scaled)[g,])
rownames(r.ft) = g

write.tsv(round(r.ft,4), "git/cluster/readRatios.tsv")

