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
# Italicizes a string.
italicize = function(x) {
  a = expression(italic(1))
  a[[1]][[2]] = x
  a
}

# Substitutes into an expression which is used as a "template",
# to simplify printing things which aren't constants (working
# around the lack of "quasiquote.")
# Args:
#   x - the expression
#   a - a list of arguments. If a name occurs as a variable in
#     the expression, its corresponding value will be substituted in.
#     (If it's not an expression, it will be coerced to one.)
# Returns: an expression, possibly modified.
expr.format = function(x, a) {

# print(x)
# print(as.character(x))
# print(class(x))

  # if it's a name, possibly substitute into it
  if (class(x) == "name") {
    x1 = as.character(x)
    if (x1 %in% names(a)) {
      a1 = a[[x1]]
      if (class(a1) == "expression")
        return(a1[[1]])
      else
        return(a1)
    }
    else {
      return(x)   # different name, so don't substitute
    }
  }

  # recursive calls, if this is an expression
  if (class(x) == "expression") {
    for(i in 1:length(x[[1]])) {
      x[[1]][[i]] = expr.format(x[[1]][[i]], a)
    }
  }

  # ... or a "call"
  if (class(x) == "call") {
    for(i in 1:length(x)) {
      x[[i]] = expr.format(x[[i]], a)
    }
  }

  x
}

# Formats a p-value, including converting smallish
# p-values to 10^{-_}.
# Args:
#   p - a p-value
#   include.equals.sign - whether to include an equals sign
#     before numbers which are not converted to 10^{-_}.
format.p = function(p, include.equals.sign = FALSE) {
  s = signif(p, 2)
  if (grepl("e", s)) {
    # XXX
    s = expression("" < 10^0)
    s[[1]][[3]][[3]][[1]] = ceiling(log10(p))
  }
  else if (include.equals.sign) {
    s = paste("=", s)
  }
  s
}

