# Various utilities for working with arrays.

# Combines arrays analogously to rbind()/cbind().
# (I'd call it "abind", but that name is already
# used, in the "abind" package.)
# Args:
#   a - a list of arrays
#   i - indices of dimensions to combine
# Returns: that array, with some dimensions combined
array.bind = function(a, i) {
  n = length(dim(a[[1]]))
  j = setdiff(c(1:n), i)    # the other dimensions

  # compute dimensions of the final array
  dims1 = rep(0, n)
  dimnames1 = as.list(rep(NA, n))
  for(k in c(1:n)) {
    if (k %in% i) {
      dims1[k] = sum(sapply(a, function(x) dim(x)[k]))
      dimnames1[[k]] = c(sapply(a, function(x) dimnames(x)[[k]]), recursive=TRUE)
    }
    else {
      s = sapply(a, function(x) dim(x)[k])
      stopifnot(all(s == s[1]))
      dims1[k] = s[1]
      # FIXME check that these are all the same?
      dimnames1[[k]] = dimnames(a[[1]])[[k]]
    }
  }
# browser()

  # shuffle these, so that dimensions to merge are first,
  # and convert to vectors
  a1 = lapply(a, function(x) as.vector(aperm(x, c(j, i))))

  # create the new array, and shuffle dims back to original order
  r1 = array(c(a1, recursive=TRUE), dim=dims1[c(j,i)])
  r = aperm(r1, order(c(j,i)))
  dimnames(r) = dimnames1
  r
}




