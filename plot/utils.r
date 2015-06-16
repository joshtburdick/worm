# Various utilities useful for plotting.

# Scales data from an interval to [0,1].
scale.interval.to.unit = function(x, interval) {
  lo = interval[1]
  hi = interval[2]

  x[ is.na(x) ] = 0
  x1 = (x - lo) / (hi - lo)
  x1[ x1 < 0 ] = 0
  x1[ x1 > 1 ] = 1
  x1
}

# Various pretty-printing utilities.
italicize = function(x) {
  a = expression(italic(1))
  a[[1]][[2]] = x
  a
}

# Formats a p-value, including converting smallish
# p-values to 10^{-_}.
format.p = function(p) {
  s = signif(p, 2)
  if (grepl("e", s)) {
    # XXX
    s = expression("" < 10^0)
    s[[1]][[3]][[3]][[1]] = ceiling(log10(p))
  }
  s
}

