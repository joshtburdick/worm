# Attempt at a particular case of unmixing.
# Efficiency isn't a priority.

source("git/unmix/ept/dirichlet.r")



# Initializes a model.
# Args:
#   m - the mixing proportions
#   x.f - the measurements in each fraction
# Returns: an initialized model
init.model = function(m, x.f) {
  num.fractions = nrow(m)
  num.cells = ncol(m)
  num.genes = nrow(x.f)
  stopifnot(ncol(x.f) == num.fractions)

  a = list(m = m, x.f = x.f,
    x = matrix(0, nrow = num.genes, ncol=num.cells))

  # add in messages from all factors
#  a$m.from.eq = 0 * a$x



  a
}






# A sample model.
sample.model = function() {
  num.frac = 2
  num.cells = 5

  # the known mixing proportions
  m = matrix(0, nrow=num.frac, ncol=num.cells)

  # the measurements (for now, just one gene)
  x.f = matrix(0, nrow=2, ncol=num.frac)

  x.f = x.f + 1
  x.f[1,1] = 10
  x.f[1,2] = 5

  m[1,1:5] = 1
  m[2,1:2] = 1

  init.model(m, x.f)
}

a = sample.model()


