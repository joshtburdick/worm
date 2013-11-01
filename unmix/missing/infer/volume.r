# Models attempting to estimate cell volume and sort matrix.

source("git/unmix/ept/approx_region_sparse.r")

# Creates a constrained linear system from assumptions
# about unmixing. This variant assumes that the amount that
# each cell is present in each fraction is known precisely
# (although it needn't be 0 or 1.)
# Args:
#   v, v.var - mean and variance of estimate of each cell's volume
#     (these should add up to about one)
#   s - estimate of how much each cell is present in each fraction
#     (between 0 and 1)
#   x.f - expression of each gene in each fraction (currently this
#     is assumed to be known exactly). Each row
# Returns: list with elements
#   A, b, b.var - constraints defining the corresponding
#     unmixing problem.
unmix.to.linear.system.1 = function(v, v.var, s, x.f) {

  # first, define the variable names
  fractions = rownames(x.f)
  genes = colnames(x.f)
  cells = names(v)
  cf = as.vector(outer(cells, fractions, 
    FUN=function(a,b) paste(a,b, sep="_")))
  gc = as.vector(outer(genes, cells,
    FUN=function(a,b) paste(a,b, sep="_")))
  v = paste(cells, "volume", sep="_")

  A.colnames = c(v, cf, gc)
  a = rep(0, length(A.colnames))
  names(a) = A.colnames

  # XXX for now, not worrying much about efficiency,
  # or making A sparse
  A = NULL
  b = NULL
  b.var = NULL

  # constraints encoding assumptions about volume
  for(i in cells) {
    a1 = a
    a1[paste(i, "volume", sep="_")] = 1
    A = rbind(A, a1)
    b = c(b, v[i])
    b.var = c(b.var, v.var[i])
  }

  # constraints on total volume in each fraction, based
  # on cell volume and the sort matrix, which we assume is
  # known exactly (here)
  for(f in fractions) {
    a1[ paste(f, "volume", sep="_") ] = -1
    for(i in cells) {
      a1 = a
      a1[ paste(i, "volume", sep="_") ] = s[f, i]

      # add this constraint
      A = rbind(A, a1)
      b = c(b, 0)
      b.var = c(b.var, 0)
    }
  }

  # constraints that expression of all genes in each cell (as a
  # fraction of the whole embryo sample) adds up to the volume of
  # each cell (again, as a fraction of the whole embryo sample)
  for(i in cells) {

    a1 = a
    a1

    for(g in genes) {


    }

    # add this constraint (again, it's exact)
    A = rbind(A, a1)
    b = c(b, 0)
    b.var = c(b.var, 0)
  }

  # constraints that expression of each gene in some fractions
  # adds up to some proportion of that fraction's total volume




  list(A = A, b = b, b.var = b.var)
}

# Creates a constrained linear system from assumptions
# about unmixing.
# Args:
#   v, v.var - mean and variance of estimate of each cell's volume
#   s, s.var - mean and variance of estimate of how much each
#     cell is present in each fraction
#   x.f, x.f.var - mean and variance of expression in each fraction
# Returns: list with elements
#   A, b, b.var - constraints defining the corresponding
#     unmixing problem.
unmix.to.linear.system = function(v, v.var, s, s.var, x.f, x.f.var) {

  # XXX for now, not worrying much about efficiency,
  # or making A sparse
  A = NULL
  b = NULL
  b.var = NULL

  # assumptions about volume


  # constraints on amount of each cell in each fraction,
  # given gene expression

  b = c(b, 0)
  b.var = c(b.var, 0)


  list(A = A, b = b, b.var = b.var)
}




