# Total coverage for each sample.

coverage.path = "git/unmix/seq/quant/rawCoverage/"
output.path = "git/unmix/seq/quant/readsPerMillion/"

# Adds up coverage for each sample, across lanes.
# Uses column names, with the first three characters
# stripped off, as the sample name.
# Args:
#   r - read counts, as a matrix with column names like, e.g.,
#     "04_cehm26"
# Returns: a similar matrix, but with columns summed by
#   sample name.
total.per.sample = function(r) {
#  sample.name = substr(colnames(r), 4, 1000)
  sample.name = sub("^0._", "", colnames(r))   # more general version of this
  r1 = t(r)
  total.counts = by(r1, sample.name, function(x) apply(x, 2, sum))
  sapply(total.counts, c)
}

# Computes reads per million mapped reads.
# FIXME this is reads per "million on the same strand", which
# is arguably incorrect.
reads.per.million = function(coverage) {
  reads.per.fraction =
    apply(coverage, 2, sum) - coverage["ribosomal_RNA",]

  reads.per.million = t( t(coverage) / (reads.per.fraction / 1e6) )

  reads.per.million
}

# Writes reads per million for some sample.
write.reads.per.million = function(coverage.file, output.file) {

  coverage.file = paste(coverage.path, coverage.file, sep="/")
  output.file = paste(output.path, output.file, sep="/")

  cat(coverage.file, "", output.file, "\n")
  r = read.table(coverage.file, sep="\t", header=TRUE,
    row.names=1, check.names=FALSE)
  rpm = reads.per.million(total.per.sample(r))
  write.table(rpm, file=output.file, sep="\t",
    row.names=TRUE, col.names=NA, quote=FALSE)
}

system(paste("mkdir -p", output.path))

if (FALSE) {
write.reads.per.million(
  "rawCoverage_20110922.tsv.gz",
  "readsPerMillion_20110922.tsv")
write.reads.per.million(
  "rawCoverage_20110922_as.tsv.gz",
  "readsPerMillion_20110922_as.tsv")

write.reads.per.million(
  "rawCoverage_Murray_050912.tsv.gz",
  "readsPerMillion_Murray_050912.tsv")
write.reads.per.million(
  "rawCoverage_Murray_050912_as.tsv.gz",
  "readsPerMillion_Murray_050912_as.tsv")

write.reads.per.million(
  "rawCoverage_Murray_52831_092812.tsv.gz",
  "readsPerMillion_092812.tsv")
write.reads.per.million(
  "rawCoverage_Murray_52831_092812_as.tsv.gz",
  "readsPerMillion_092812_as.tsv")
}

write.reads.per.million(
  "rawCoverage_embryo_ts.tsv.gz",
  "readsPerMillion_embryo_ts.tsv")


