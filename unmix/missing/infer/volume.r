# Models attempting to estimate expression and cell volume
# (and possibly the sort matrix.)

source("git/unmix/ept/approx_region_sparse.r")

# Creates a constrained linear system from assumptions
# about unmixing. This variant assumes that the amount that
# each cell is present in each fraction is known precisely
# (although it needn't be 0 or 1.)
# Args:
#   v, v.var - mean and variance of estimate of each cell's volume
#     (these should add up to about one)
#   s - estimate of how much each cell is present in each fraction
#     (betwen 0 and 1, I think)
#   x.f - expression of each gene in each fraction (currently this
#     is assumed to be known exactly). Each row corresponds to one
#     gene, and each column corresponds to one row of s.
# Returns: list with elements
#   A, b, b.var - constraints defining the corresponding
#     unmixing problem.
unmix.to.linear.system.1 = function(v, v.var, s, x.f) {

  # first, define the variable names
  # ??? do I need these?
  fractions = rownames(x.f)
  genes = colnames(x.f)
  cells = names(v)
  cf = as.vector(outer(cells, fractions, 
    FUN=function(a,b) paste(a,b, sep="_")))
  gc = as.vector(outer(genes, cells,
    FUN=function(a,b) paste(a,b, sep="_")))
  v = paste(cells, "volume", sep="_")

  # tables storing constraints (in sparse form)
  A1 = NULL
  b1 = NULL

  # constraints encoding assumptions about volume of each cell
  for(i in cells) {

    A1 = rbind(A1,
      data.frame(i = paste(cells, "volume_assumption", sep="_"),
        j = paste(cells, "volume", sep="_"),
        a = 1, stringsAsFactors=FALSE))

    b1 = rbind(b1,
      data.frame(i = cells,
        b = v[cells],
        b.var = v.var[cells], stringsAsFactors=FALSE))
  }

  # constraints on total volume in each fraction, based on cell volume
  # and the sort matrix, which we assume is known exactly (here).
  # possibly not needed
  for(f in fractions) {

    A1 = rbind(A1,
      data.frame(i = paste(f, "volume", sep="_"),
        j = paste(f, cells, sep="_"),
        a = s[f, cells], stringsAsFactors=FALSE))
    A1 = rbind(A1,
      data.frame(i = paste(f, "volume", sep="_"),
        j = paste(f, "volume", sep="_"),
        a = -1, stringsAsFactors=FALSE))

    b1 = rbind(b1,
      data.frame(i = paste(f, "volume", sep="_",
        b = 0,
        b.var = 0, stringsAsFactors=FALSE))
  }

  # constraints that total expression of all genes in each cell
  # adds up to the volume of that cell (both as a fraction of
  # the whole embryo sample)
  # also possibly not needed
  for(i in cells) {

    A1 = rbind(A1,
      data.frame(i = ,
        j = paste(


  }

  # constraints that expression of each gene in each fractions
  # adds up to some proportion of that fraction's total volume
  # (currently these are assumed to be exact)
  for(f in fractions)
    for(g in genes) {

      A1 = rbind(A1,
        data.frame(i = paste(g, f, sep="_"),
          j = paste(g, cells, sep="_"),
          a = s[f, cells], stringsAsFactors=FALSE))

      b1 = rbind(b1, 
        data.frame(i = ,
          b = 0,
          b.var = 0, stringsAsFactors = FALSE))

    }

  list(A1 = A1, b1 = b1)
#  list(A = A, b = b, b.var = b.var)
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




