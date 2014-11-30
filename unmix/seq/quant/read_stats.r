# Statistics about total number of reads mapped.

# total reads mapped
mapped.reads = as.matrix(read.table(gzfile(
  "R/unmix/sort_paper/seq/quant/rawCoverage.tsv.gz"),
    sep="\t", header=TRUE, row.names=1))
mapped.reads.per.sample = apply(mapped.reads, 2, sum)
s = names(mapped.reads.per.sample)
s = sub("\\.tsv\\.gz", "", s)
s = substr(s, 5, 1000)
total.mapped = c(by(mapped.reads.per.sample, s, sum))

read.stats = data.frame(total.mapped, stringsAsFactors=FALSE)
write.table(read.stats, file="git/unmix/seq/quant/read_stats.tsv",
  sep="\t", row.names=TRUE, col.names=NA)

