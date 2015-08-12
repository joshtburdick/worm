# Various hopefully-useful utilities.

# Reads a .tsv file (should read files produced
# by write.tsv(), below.)
read.tsv = function(file)
  read.table(file=file, sep="\t", header=TRUE, row.names=1,
    check.names=FALSE, as.is=TRUE)

# Writes a .tsv file in a format that Excel can read.
write.tsv = function(x, file)
  write.table(x, file=file, sep="\t", na="",
    row.names=TRUE, col.names=NA, quote=FALSE)

# useful for printing progress
backspace.string = paste(rep("\b", 60), collapse="")

# a similar function
write.status = function(s) {
  max.length = 70
  s = substr(s, 1, max.length)

  bs = paste(rep("\b", max.length), collapse="")
  pad = paste(rep(" ", max.length - nchar(s)), collapse="")

  cat(paste0(bs, s, pad))
}

# Constructs an array, with dimensions implied by dimnames.
# Args:
#   dimnames - a list giving the dimnames of the array
# Returns: an array, with those dimnames.
array.from.dimnames = function(dimnames) {
  array(dim = sapply(dimnames, length),
    dimnames = dimnames)
}

# Simple "named capture" function.
# Args:
#   s - a character vector
#   pattern - a regexp
# Returns: substring matching that pattern.
regexp.capture = function(s, pattern)
  regmatches(s, gregexpr(pattern, s))[[1]]

