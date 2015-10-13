# Simulates expression with a gamma distribution.

# Converts from "number of lineages on" to a random
# gamma-distributed expression value.
# Args:
#   shape.m, scale.m
#   a - matrix of how many lineages are on in a given cell
# Returns: matrix of simulated expression data
sim.expr.gamma.1 = function(shape.m, scale.m, a) {
  x = rgamma(length(a), shape=shape.m*a+1, scale=scale.m*a+1)

  x = array(x, dim = dim(a))
  dimnames(x) = dimnames(a)

  rownames(x) = paste(rownames(x), shape.m, scale.m, sep="_")
#  x = x / apply(x, 1, sum)
  x
}

# Given a matrix of on-off cells, creates a data matrix with
# gamma-distributed expression based on that matrix.
sim.expr.gamma = function(a) {
  set.seed(42)
  list(on.off = a,
    x = rbind(
      sim.expr.gamma.1(2, 0, a),
      sim.expr.gamma.1(2, 2, a),
      sim.expr.gamma.1(2, 4, a),
      sim.expr.gamma.1(2, 9, a),
      sim.expr.gamma.1(4, 0, a),
      sim.expr.gamma.1(4, 2, a),
      sim.expr.gamma.1(4, 4, a),
      sim.expr.gamma.1(4, 9, a),
      sim.expr.gamma.1(9, 0, a),
      sim.expr.gamma.1(9, 2, a),
      sim.expr.gamma.1(9, 4, a),
      sim.expr.gamma.1(9, 9, a)))
}

