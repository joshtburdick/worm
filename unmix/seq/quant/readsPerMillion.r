# Total coverage for each sample.

# Adds up coverage for each sample, across pools.
# Args:
#   r - read counts, as a data frame with column names like, e.g.,
#     "X04_cehm26.tsv.gz"
# XXX this makes many assumptions about the number of lanes, sample names, etc.
# Returns: total read counts for each sample (with slightly
#   modified column names)
# FIXME probably should use this in compute.reads.per.million().
total.per.pool = function(r) {
  fraction.name = gsub("^X0._", "", gsub(".tsv.gz", "", colnames(r)))

  total.coverage = matrix(0, nrow=nrow(r), ncol=24)
  rownames(total.coverage) = rownames(r)
  colnames(total.coverage) = fraction.name[c(1:12, 37:48)]

  for(i in 1:ncol(r))
    total.coverage[,fraction.name[i] ] = total.coverage[,fraction.name[i] ] + r[,i]

  replace.m.with.hyphen = function(a) {
    ch = substr(a, 4, 4)
    ch[ch=="m"] = "-"
    paste(substr(a, 1, 3), ch, substr(a, 5, 1000), sep="")
  }

  colnames(total.coverage) = replace.m.with.hyphen(colnames(total.coverage))
  total.coverage
}

# Hack to average samples with the same name.
average.by.sample.name = function(r)
  cbind(
  "all" = (r[,"cnd-1_singlets"] + r[,"pha-4_singlets"]) / 2,
  "ceh-26" = r[,"ceh-26"],
  "ceh-27" = r[,"ceh-27"],
  "ceh-36" = r[,"ceh-36p"],
  "ceh-36_minus" = r[,"ceh-36m"],
  "ceh-6" = r[,"ceh-6"],
  "cnd-1" = (r[,"cnd-1_12_14"] + r[,"cnd-1p1_4"] + r[,"cnd-1p8_19"]) / 3,
  "cnd-1_minus" = r[,"cnd-1m"],
  "F21D5.9" = r[,"F21D5.9"],
  "hlh-16" = r[,"hlh-16"],
  "irx-1" = r[,"irx-1"],
  "mir-57" = r[,"mir-57"],
  "mls-2" = r[,"mls-2"],
  "pal-1" = r[,"pal-1"],
  "pha-4" = (r[,"pha-4p12_9"] + r[,"pha-4p9_1"]) / 2,
  "pha-4_minus" = r[,"pha-4m"],
  "ttx-3" = r[,"ttx-3"],
  "unc-130" = r[,"unc-130"])

# Computes reads per million mapped reads.
# FIXME this is reads per "million on the same strand", which
# is arguably incorrect.
get.reads.per.million = function(coverage.file) {

  r = read.table(gzfile(coverage.file), sep="\t", header=TRUE, row.names=1)

  fraction.name = gsub("^X0._", "", gsub(".tsv.gz", "", colnames(r)))
  total.coverage = matrix(0, nrow=nrow(r), ncol=24)
  rownames(total.coverage) = rownames(r)
  colnames(total.coverage) = fraction.name[c(1:12, 37:48)]

  for(i in 1:ncol(r))
    total.coverage[,fraction.name[i] ] = total.coverage[,fraction.name[i] ] + r[,i]

  reads.per.fraction = apply(total.coverage, 2, sum) - total.coverage["ribosomal_RNA",]
  reads.per.million = t( t(total.coverage) / (reads.per.fraction / 1e6) )

  replace.m.with.hyphen = function(a) {
    ch = substr(a, 4, 4)
    ch[ch=="m"] = "-"
    paste(substr(a, 1, 3), ch, substr(a, 5, 1000), sep="")
  }
  colnames(reads.per.million) = replace.m.with.hyphen(colnames(reads.per.million))

  reads.per.million
}

reads.per.million = get.reads.per.million("git/unmix/seq/quant/rawCoverage.tsv.gz")

write.table(round(average.by.sample.name(reads.per.million), 3),
  file=gzfile("git/unmix/seq/quant/readsPerMillion.tsv.gz"), sep="\t", col.names=NA)


