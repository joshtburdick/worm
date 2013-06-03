# Definition of simulated expression when a cell is "on".

# Converts from "number of lineages on" to a random expression value.
# Args:
#   a - matrix of how many lineages are on in a given cell
# Returns: matrix of simulated expression data,
sim.expr = function(a) {
  x = rnorm(length(a), mean=10*a+1, sd=sqrt(a+1))  # was 50*a+1
  x[ x<0 ] = 0
  x = array(x, dim = dim(a))
  dimnames(x) = dimnames(a)
#  x = x / apply(x, 1, sum)
  x
}


