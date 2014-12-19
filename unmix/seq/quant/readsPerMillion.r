# Total coverage for each sample.

coverage.path = "git/unmix/seq/quant/rawCoverage/WS220_20140111/"
output.path = "git/unmix/seq/quant/readsPerMillion/WS220_20140111/"

# Strip off lane number, if it's present, from some names.
remove.lane.name = function(a) sub("^0._", "", a)

# Adds up coverage for each gene in each sample, across lanes.
# Uses column names, with the first three characters
# potentially stripped off, as the sample name.
# Args:
#   r - read counts, as a matrix with column names like, e.g.,
#     "04_cehm26"
# Returns: a similar matrix, but with columns summed by
#   sample name.
gene.total.per.sample = function(r) {
  # strip off lane number, if it's present
  sample.name = remove.lane.name(colnames(r))
  r1 = t(r)
  total.counts = by(r1, sample.name, function(x) apply(x, 2, sum))
  sapply(total.counts, c)
}

# Computes total reads per sample (based on the counts written by
# "get_coverage.pl", which will be somewhat larger than the sum
# over all genes.)
# Args:
#   rt - the read totals
# Returns: 
total.per.sample = function(rt, r) {
  c(by(rt[,1], remove.lane.name(rownames(rt)), sum))
}

# Computes reads per million mapped reads.
# Args:
#   coverage - the read coverage for each gene
#   total.reads - the total number of mapped reads for each gene
# Returns: list with
#   rpm - reads per million of each gene
#   rt - the adjusted read totals for each sample
reads.per.million = function(coverage, total.reads) {

  # the same, with rRNA and mitochondrial DNA subtracted out
  genes.to.ignore = grep("^(rrn-|MTCE.)", rownames(coverage), value=TRUE) 
  print(genes.to.ignore)
  reads.to.ignore = apply(coverage[genes.to.ignore,], 2, sum)
  reads.per.fraction =
    total.reads[ colnames(coverage) ] - reads.to.ignore[ colnames(coverage) ]

  reads.per.million = t( t(coverage) / (reads.per.fraction / 1e6) )

  g = setdiff(rownames(coverage), genes.to.ignore)
  list(rpm = reads.per.million[ g, ], rt = reads.per.fraction)
}

# Writes reads per million for some sample.
write.reads.per.million = function(coverage.file, output.basename) {

  coverage.file = paste(coverage.path, coverage.file, sep="/")

  cat(coverage.file, "\n")
  r = read.table(coverage.file, sep="\t", header=TRUE,
    row.names=1, check.names=FALSE)
  rt = read.table(
    sub(".tsv.gz", ".totalReads.tsv", coverage.file),
    header=TRUE, row.names=1, check.names=FALSE)

  a = reads.per.million(gene.total.per.sample(r), total.per.sample(rt))
  write.table(a$rpm, file=paste0(output.path, "/", output.basename, ".tsv"),
    sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
  write.table(a$rt, file=paste0(output.path, "/", output.basename, "_read_totals.tsv"),
    sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
}

system(paste("mkdir -p", output.path))

if (TRUE) {

write.reads.per.million(
  "Murray050912.tsv.gz",
  "050912")
write.reads.per.million(
  "Murray050912_as.tsv.gz",
  "050912_as")

write.reads.per.million(
  "Murray092812.tsv.gz",
  "092812")
write.reads.per.million(
  "Murray092812_as.tsv.gz",
  "092812_as")
}

if (TRUE) {
write.reads.per.million(
  "20110922.tsv.gz",
  "20110922")
write.reads.per.million(
  "20110922_as.tsv.gz",
  "20110922_as")
}

write.reads.per.million(
  "embryo_ts.tsv.gz",
  "embryo_ts")

