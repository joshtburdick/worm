# Another attempt at estimating expression, and the sort matrix,
# simultaneously.

source("git/utils.r")
# source("git/unmix/ept/gamma.r")
source("git/unmix/ml/pos_linear_solve.r")
source("git/unmix/ml/unmixLikelihood.r")

# utilities, slow-ish but possibly useful
gamma.s2mv = function(a) list(m = a$a / a$b, v = a$a / (a$b^2))
gamma.mv2s = function(a) {
  b = a$m / a$v
  list(a = a$m * b, b = b)
}
mv2moment = function(a) list(x1 = a$m, x2 = a$m^2 + a$v)
moment2mv = function(a) list(m = a$x1, v = a$x2 - (a$x1^2))

# Does one step of "updating", pulling a point towards the nearest
# point which satisfies the constraints (but stopping where
# any entry would go to 0.)
# Args:
#   A, B - these give the linear constraint
#   X - this gives the current solution
# Returns: an updated X, which is closer to a solution of
#   AX = B, with X >= 0.
pos.linear.solve.1 = function(A, B, X) {
  X1 = lin.constraint(A, B)(X)
  X = move.pos(X, X1) + 1e-20
  X
}

# Normalizes the rows of a matrix.
norm.rows = function(a) a / apply(a, 1, sum)

# "Updates" a gamma distribution with one observation.
# Args:
#   p - a gamma prior
#   x - a data point
# Returns: mean and variance of posterior.
gamma.update = function(p, x) {
  a = mv2moment(gamma.s2mv(p))
  b = list(x1 = (a$x1 + x) / 2, x2 = (a$x2 + x^2) / 2)
  moment2mv(b)$m
}

# Estimates expression, and the sort matrix.
# Args:
#   a.prior - prior on the sort matrix
#   x.prior - prior on expression
#   b - read data, as proportions (each row should sum to 1)
#   max.iters - maximum number of iterations to do
#   save.hist - whether to save a history of estimates
# Returns: list with elements:
#   a - estimated sort matrix
#   x - estimated expression
#   update.stats - likelihood stats, and how much a and x changed
unmix.expr.and.sort.matrix.1 =
    function(a.prior, x.prior, b, max.iters=50, save.hist=TRUE) {

  # these track several measures of how good the solution is
  update.stats = NULL
  h = list()

  # initialize estimates to the means of the priors
  # (the constraint will admittedly be way off, initially)
  a = gamma.s2mv(a.prior)$m
  x = gamma.s2mv(x.prior)$m

  for(iter in 1:max.iters) {

    # update expression
    x1 = pos.linear.solve.1(a, b, x)
    x2 = gamma.update(x.prior, t(norm.rows(t(x1))))$m

    # update sort matrix
    a1 = pos.linear.solve(t(x2), t(b), a)
    a2 = gamma.update(a.prior, norm.rows(a1))$m

    # save updates
    x = x2
    a = a2

    # collect stats about how "good" the solution is
    update.stats.1 = unmix.ll.and.fit(
      list(a.prior = a.prior, x.prior = x.prior, b = b),
      list(a = a, x = x))
write.status(paste(update.stats.1, collapse=" "))
    update.stats = rbind(update.stats, update.stats.1)
    if (save.hist) {
      h[[ iter ]] = list(a = a, x = x)
    }
  }

  list(a = a, x = x, update.stats = update.stats, history = h)
}

