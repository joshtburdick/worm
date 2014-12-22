# Attempt at a more general library for EP.

# A model is represented by a list with at least:
#   vars - a named list of variables (the names should
#     be unique)
#   factors - a list of factors

# Makes a model.
# Args:
#   vars - list of variables (this is just a list of
#     exponential-family variables)
#   factors - list of factors (see make.factors())
make.model = function(vars, factors) {
  # XXX this admittedly looks pretty silly
  r = list(vars = vars, factors = factors)
#  class(x) = "epModel"
  r
}

# XXX this is deprecated
# A variable is a list with at least these elements:
#   b - the total messages to this variable (as a matrix)
#   observed - TRUE iff this variable is observed
#   log.evidence - function which computes log-evidence
# Makes a variable.
make.var.old = function(b, log.evidence) {
  r = list(b = b, observed = FALSE, log.evidence = log.evidence)
#  attr(r, "class") = "ep.var"
  r
}

# Adds a scaled amount of some variable to another variable.
# Possibly not used.
# Args:
#   x, y - lists of numeric variables (y should only have
#     a subset of x's variables)
#   a - scalar
# Returns: x + ay
add.scaled = function(x, y, a) {
  for(i in names(y))
    x[[i]] = x[[i]] + a * y[[i]]
  x
}

# Adds in a message.
message.add = function(x, y) {
  for(i in names(y))
    x[[i]] = x[[i]] + y[[i]]
  x
}

# Subtracts out a message.
message.subtract = function(x, y) {
  for(i in names(y))
    y[[i]] = x[[i]] - y[[i]]
  y
}

# A factor is a list with:
#   m - list of messages from this factor
#   f - list of methods, including:
#     update - function which, given messages to this factor,
#       computes messages from this factor
#     logEvidence - function which computes log-evidence for this factor

# Constructs a factor with some set of methods.
# Args:
#   x - the variables connected to this factor
#   f - list containing the "update" and "logEvidence" functions
#     for this factor
make.factor = function(x, f) {
  # the messages start at the belief, but zeroed out
  z = lapply(x, function(a) 0 * a$b )
  names(z) = names(x)

  a = list(m = z, f = f)
#    update = f$update, logEvidence = f$logEvidence)
#  class(a) = "ep.factor"
  a
}

# An equality factor.
eq.f = function(a) {
  a1 = a[[1]]
  a2 = a[[2]]
  
  list(update = list(a2, a1),
    logEvidence = logPartition(a1+a2) - (logPartition(a1) + logPartition(a2)))
}

# Factor which constrains a variable (almost the same thing.)
# ??? should this care about the message passed to it?
constrain.f = function(x) function(a) {
  a1 = a[[1]]

  list(update = list(x),
    logEvidence = log.partition(x+a1) - (log.partition(x) + log.partition(a1)))
}

# Updates a model once, using "parallel" updating.
# Currently, no convergence testing.
#   m - a model, as a list with elements "vars" and "factors"
# Returns: m, after one round of parallel updates
update.ep = function(m) {

  f1 = list()  

  for(i in 1:length(m$factors)) {

    # message to that factor
    m.to = message.subtract(m$vars, m$factors[[i]]$m)

    # compute updated messages from this factor
    m.from = m$factors[[i]]$f$update(m.to)

    f1[[i]] = list(m = m.from, f = m$factors[[i]]$f)
  }

  # clear the variables
  v1 = lapply(m$vars, function(v) 0 * v)

  # send messages from each factor
  for(i in 1:length(f1)) {
    v1 = message.add(v1, f1[[i]]$m)
  }

  list(vars = v1, factors = f1)
}

