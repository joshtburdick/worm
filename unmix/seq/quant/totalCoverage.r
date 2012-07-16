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

# Computes reads per million mapped reads.
# FIXME this is reads per "million on the same strand", which
# is arguably incorrect.
compute.reads.per.million = function(coverage.file, output.file) {

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

  write.table(round(reads.per.million, 3), file=gzfile(output.file), sep="\t", col.names=NA)
}

compute.reads.per.million("git/unmix/seq/quant/rawCoverage.tsv.gz",
  "git/unmix/seq/quant/readsPerMillion.tsv.gz")


