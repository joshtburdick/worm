# Summary statistics about the number of reads mapped.

source("git/utils.r")

r = read.tsv("git/unmix/seq/quant/read_stats.tsv")$mapped.reads.facs

cat("total reads =", sum(r), "\n")
cat("mean reads per sample =", mean(r), "\n")
cat("range of reads per sample =", min(r), "to", max(r), "\n")
cat("median reads per sample =", median(r), "\n")


