# Adds a (currently ignored) "score" column to
# upstream ChIP counts. (If the "score" column is
# present, then this simply returns it unmodified.)

# Possibly not used.

add.score = function(counts) {
  if ("score" %in% names(dimnames(counts)))
    return(counts)

  stopifnot(dim(counts) == 3)

  dim1 = c(dim(counts), 1)
  dimnames1 = c(dimnames(motif.count), list(score=c("0")))

  array(counts, dim=dim1, dimnames=dimnames1)
}

