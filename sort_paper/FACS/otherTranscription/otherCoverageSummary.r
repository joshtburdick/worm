# Summary of coverage outside of annotated genes.

source("git/utils.r")

result.base = "git/sort_paper/FACS/otherTranscription/otherCoverage/"

# Histogram data by log10 buckets.
# Args:
#   x - the data to histogram
#   num.buckets - how many buckets to include per factor-of-10
#     increase in x
#   max.x - where to clip the numbers.
# Returns: a vector of counts of log10 numbers.
log10.histogram = function(x, num.buckets = 4, max.x = 1e6) {

  # trim data
  x[ x < 1 ] = 1
  x[ x >= max.x ] = max.x

  # create the buckets
  buckets = c(0:floor(num.buckets * log10(max.x)))
  r = rep(0, length(buckets))
  names(r) = buckets

  # bucketize, and add up counts
  counts = c(table(floor(num.buckets * log10(x))))
  r[ names(counts) ] = counts
  r
}

# Adds up total coverage in a given .bed file.
# Args:
#   f - a filename
# Returns: histogram of number of reads, and size of group
# of reads, for that file.
coverage.histogram = function(f) {
  r = read.table(f, sep="\t", header=FALSE, as.is=TRUE)
  colnames(r) = c("chr", "a", "b", "numReads", "strand")
  r$size = r$b - r$a
  list(reads = log10.histogram(r$numReads),
    size = log10.histogram(r$size))
}

r = NULL

# omit HS, N2, and RNAi samples
files1 = list.files(paste0(result.base, "/outside_exons/"))
files2 = grep("^(HS|N2|ges1|lit1|pop1)", files1, value=TRUE, invert=TRUE)


for (f in files2) {
  a = sub(".bed", "", f)
  write.status(a)
  exons = coverage.histogram(paste0(result.base, "/outside_exons/", f))
  genes = coverage.histogram(paste0(result.base, "/outside_genes/", f))
  r = rbind(r, c(sample = a, stat = "read count outside exons", exons$reads))
  r = rbind(r, c(sample = a, stat = "read count outside genes", genes$reads))
  r = rbind(r, c(sample = a, stat = "group size outside exons", exons$size))
  r = rbind(r, c(sample = a, stat = "group size outside genes", genes$size))
}

write.tsv(r, "git/sort_paper/FACS/otherTranscription/otherCoverageSummary.tsv")

