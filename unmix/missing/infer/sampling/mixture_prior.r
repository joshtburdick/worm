# Informative prior for the mixture proportions.

library(corpcor)

source("git/utils.r")

# Unnormalized log-Dirichlet density.
ld.dirichlet = function(alpha, x) {
  sum((alpha - 1) * log(x))
}

# Log-likelihood for the sorting model.
# Args:
#   p - the prior, as a list with Dirichlet parameters:
#     V.alpha - the volume
#     S.alpha - how each cell was sorted. An array indexed by:
#       fraction, cell, and where sorted ("plus" or "minus")
#   x - the corresponding values, as a list with
#     elements V and S
sort.dirichlet.ll = function(p, x) {
  ld.dirichlet(p$V.alpha, x$V) +
    ld.dirichlet(p$S.alpha[,,"plus"], x$S) +
    ld.dirichlet(p$S.alpha[,,"minus"], 1-x$S) 
}

# Normalizes a sort matrix to the form actually used
# in unmixing.
sort.matrix.normalize = function(x) {
  # compute "proportion of embryo in each fraction"
  m.plus = t( x$V * t(x$S) )
  m.minus = t( x$V * (1 - t(x$S)))

  rownames(m.plus) =
    paste(rownames(m.plus), " (+)", sep="")
  rownames(m.minus) =
    paste(rownames(m.minus), " (-)", sep="")
  m = rbind(m.plus, m.minus)
#  colnames(m) = FIXME

  # FIXME: add in doubly-sorted fractions

  # normalize the rows
  m = m / apply(m, 1, sum)
  m
}

# A prior on the sort matrix. This version of it assumes
# uncertainty in sorting, and proportion of cells missing,
# is consistent across the lineage; however, the underlying
# model allows being more specific.
# Args:
#   V.alpha - the cell volume Dirichlet parameters
#   m - the sort matrix (with rows unnormalized); entries
#     should be between 0 and 1
#   concentration - concentration for the sorting beta
#     (for now, this is constant across the sort matrix)
# Returns: a list of parameters (see sort.dirichlet.ll):
#   V.alpha - the volume
#   S.alpha - how each cell was sorted, as an array indexed by:
#     plus, minus (and later, perhaps, missing)
sort.prior.dirichlet.1 = function(V.alpha, m, concentration) {
  names(V.alpha) = colnames(m)

  S.alpha = array(dim=c(nrow(m), ncol(m), 2))
  dimnames(S.alpha) = list(fraction=rownames(m),
    cell=colnames(m),
    sort=c("plus", "minus"))   # ??? add "missing"?

  S.alpha[,,"plus"] = concentration * m
  S.alpha[,,"minus"] = concentration * (1-m)

  list(V.alpha = V.alpha, S.alpha = S.alpha)
}

# Picks an initial sample from the prior.
sort.init = function(p) {
  V = rgamma(length(p$V.alpha), shape=p$V.alpha, scale=1)
  V = V / sum(V)
  names(V) = names(p$V.alpha)  

  S = matrix(rbeta(length(p$S.alpha[,,"plus"]),
    shape1 = p$S.alpha[,,"plus"], shape2 = p$S.alpha[,,"minus"]),
    nrow = nrow(p$S.alpha[,,"plus"]))
  rownames(S) = dimnames(p$S.alpha)[[1]]
  colnames(S) = dimnames(p$S.alpha)[[2]]

  list(V = V, S = S)
}

# Converts from V and S to the sort matrix.
# Args:
#   V, S - cell volumes, and amount each cell is sorted
# Returns: the corresponding sort matrix
to.sort.matrix = function(V, S) {
  m = t( V * t(S) )
  m / apply(m, 1, sum)
}

# Linear constraint, for a matrix of numbers.
# Args:
#   A, B - these define the constraint
#   m - the current mean
# Returns: the point nearest to m which satisfies
#   the constraint.
lin.constraint.matrix = function(A, B) {

  P = pseudoinverse( A %*% t(A) )

  # we return a function which computes this
  function(M) {
#    browser()
    M - t(A) %*% P %*% (A %*% M - B)
  }
}

# Goes "as far as possible" in some direction, without
# any element of a vector going negative.
move.bounded.matrix = function(a, b) {
  a = t(a)
  b = t(b)

  # the amount we would go in each direction
  delta = b - a

  s = - a / delta
  s[ is.na(s) ] = 0

  # If any element of s is negative, it doesn't matter, in
  # terms of hitting boundaries. However, if any element is
  # positive, we can only add that "amount" of delta.
  # was: m = min(1, s[ s > 0 ])
  s[ s < 0 ] = 1
  s[ s > 1 ] = 1
  m = apply(s, 1, min)

  t( a + m * delta )
}

# Picks an initial sample of the expression
# (either a point satisfying the constraints, or
# a point as close as possible to that.)
# Args:
#   A, B - these define the linear constraints
# Returns: list with elements
#   X - expression matrix such that
#     AX = B, X >= 0, each row of X sums to 1
#     (or as near as possible to that)
#   del - the difference in the last update
expr.init = function(A, B) {

  lc = lin.constraint.matrix(A, B)

  # initially, x is the pseudoinverse solution
  x = lc(matrix(0, nrow=ncol(A), ncol=ncol(B)))
  del = Inf

  # iteratively refine the solution
  for(iter in 1:3) {
    x1 = x

    # enforce the linear constraint
    x = lc(x)

    # and positivity (FIXME: use move.bounded() ? )
    # x[ x < 0 ] = 0
    x = move.bounded.matrix(x1, x)

    # and rows summing to 1
    x = x / apply(x, 1, sum)

    del = sum(x1 - x)
    write.status(paste(iter, del))
  }

  list(x = x, del = del)
}

# Perturbs the sort matrix somewhat. 
sort.matrix.jump.1 = function(x) {
print(range(x$S, na.rm=TRUE))
  x$S[ is.na(x$S) ] = 0
  V = x$V * rgamma(length(x$V), shape=50, scale=0.02)
  V = V / sum(V)
  S1 = rbeta(nrow(x$S) * ncol(x$S), 1e-3 + 50 * x$S, 1e-3 + 50 * (1-x$S))

  list(V = V, S =
    matrix(S1, nrow=nrow(x$S), dimnames=dimnames(x$S)))
}

# Tries to find an initial solution by optimization.
# Args:
#   m - the sort matrix parameters (as V and S)
#   b - the expression measurements
#   x - the initial expressionn
find.init = function(m, b, num.iters = 10) {
  err.best = Inf
  r = NULL
  x = NULL

  for(i in 1:num.iters) {

    m1 = sort.matrix.jump.1(m)
    m1.n = sort.matrix.normalize(m)
    sort.fractions = intersect(rownames(m1.n), rownames(b))
    m1.n = m1.n[ sort.fractions , ]

    r = expr.init(m1.n, b)

    # compute error in the constraint
    z = m1.n %*% r$x - b
    err = sum(z^2)
    cat(paste("\n", i, err, "\n"))

    # possibly update "best"
    if (err < err.best) {
      m = m1
      x = r$x
      err.best = err
    }
  }

  list(m = m, x = x, err.best = err.best)
}

