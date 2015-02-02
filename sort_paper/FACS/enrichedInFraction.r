# Computes which genes were enriched in which fraction.

source("git/utils.r")

r = read.tsv("git/cluster/readRatios.tsv")

r.sort.only = r[,c(1:23)]

r.sort.only.averaged = as.matrix(cbind(r.sort.only[,c(7:21)],
  "cnd-1" = apply(r.sort.only[,c(1:3)], 1, mean),
  "pha-4" = apply(r.sort.only[,c(4:6)], 1, mean),
  "singlets" = apply(r.sort.only[,c(22:23)], 1, mean)))

r.sort.only.averaged[ is.na(r.sort.only.averaged) ] = 0

# Given a data table, constructs a list of things
# enriched and depleted.
# Args:
#   r - a data table
#   cutoff - the cutoff to use
# Returns: a list of boolean vectors, one per column of r.
get.enriched.and.depleted = function(r, cutoff) {
  a = NULL
  for(s in colnames(r)) {
    # XXX slightly confusing-looking code
    a[[ paste(s, "enriched") ]] = (r[,s] >= cutoff)[ r[,s] > -cutoff ]
    a[[ paste(s, "depleted") ]] = (r[,s] <= -cutoff)[ r[,s] < cutoff ]
  }
  a
}

facs.enriched.depleted = get.enriched.and.depleted(r.sort.only.averaged, 2)


