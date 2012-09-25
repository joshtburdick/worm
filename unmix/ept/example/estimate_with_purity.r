# Attempt at a particular case of unmixing.
# Efficiency isn't a priority.

source("git/unmix/ept/dirichlet.r")

num.frac = 2
num.cells = 5

# a sample model; completely a toy example
a = list(

  # the measurements
  f = matrix(0, nrow=2, ncol=num.frac),

  # the known mixing proportions
  m = matrix(0, nrow=num.frac, ncol=num.cells),

  # the expression estimates
  x = matrix(0, nrow=2, ncol=num.cells)
)

a$f = a$f + 1
a$f[1,1] = 10
a$f[1,2] = 5

a$m[1,1:5] = 1
a$m[2,1:2] = 1

# add in messages from all factors
a$m.from = 0 * a$x








