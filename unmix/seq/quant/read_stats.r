# Statistics about total number of reads mapped.

# Total reads mapped for one set of samples.
# Args:
#   raw.coverage.file - file containing the raw coverage
#   trim.names - if TRUE, remove the initial \
# Returns: vector of read counts
count.mapped = function(raw.coverage.file, trim.names=TRUE) {
  coverage.dir = "git/unmix/seq/quant/rawCoverage/WS220_20140111/"
  mapped.reads = as.matrix(read.table(gzfile(
    paste0(coverage.dir, raw.coverage.file)),
      sep="\t", header=TRUE, row.names=1))
  mapped.reads.per.sample = apply(mapped.reads, 2, sum)
  s = names(mapped.reads.per.sample)
  s = sub("\\.tsv\\.gz", "", s)

  # XXX this removes the lane number
  if (trim.names) {
    s = substr(s, 5, 1000)
  }
  total.mapped = c(by(mapped.reads.per.sample, s, sum))
  total.mapped = total.mapped[ sort(names(total.mapped)) ]
  total.mapped
}

total.mapped = c(
  count.mapped("20110922.tsv.gz", trim.names=FALSE) +
    count.mapped("20110922_as.tsv.gz", trim.names=FALSE),
  count.mapped("Murray050912.tsv.gz") +
    count.mapped("Murray050912_as.tsv.gz"),
  count.mapped("Murray092812.tsv.gz") +
    count.mapped("Murray092812_as.tsv.gz"))

# subset of read stats for just the FACS-sorted experiments
mapped.reads.facs = total.mapped[ grep("^(ges1_|lit1_|pop1_|HS|N2_)", names(total.mapped), value=TRUE, invert=TRUE) ]

# XXX for now, just including the reads from the FACS sorting
read.stats = data.frame(mapped.reads.facs, stringsAsFactors=FALSE)
# write.table(read.stats, file="git/unmix/seq/quant/read_stats.tsv",
#   sep="\t", row.names=TRUE, col.names=NA)

