# Utilities using biomaRt.

library("biomaRt")

# wb.gene = useMart("WS220", "wormbase_gene")

# Generalized "conversion" function.
# Args:
#   mart - a biomaRt object
#   from - attribute to filter on
#   to - attribute to return
# Returns: function which converts names.
mart.convert = function(mart, from, to) function(x) {
  r = getBM(attributes = c(from, to),
    filters = c(from),   
    values = list(x),
    mart = mart)
  r
  r1 = r[ match(x, r[,1]) , 2 ]
  r1
}

