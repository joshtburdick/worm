# Computes enrichment of things in e.g. enriched fractions versus depleted,
# comparing, e.g., things enriched versus things depleted.

source("git/sort_paper/FACS/enrichedInFraction.r")

# Sets of genes enriched, versus depleted.
# Args:
#   r - a data table
#   cutoff - the cutoff to use
# Returns: a list of boolean vectors, one per column of r.
enriched.vs.depleted.1 = function(r, cutoff) {
  a = NULL
  for(s in colnames(r)) {
    a[[ s ]] =
      (r[,s] >= cutoff)[ abs(r[,s]) > cutoff ]
  }
  a
}

enriched.vs.opposite.gene.sets = {
  plus = enriched.vs.depleted.1(r.sort.only.averaged, 1)
  names(plus) = paste(names(plus), "enriched")
  minus =  enriched.vs.depleted.1(-r.sort.only.averaged, 1)
  names(minus) = paste(names(minus), "depleted")
  c(plus, minus)
}



