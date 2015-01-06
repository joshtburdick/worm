# Attempt at a more general library for EP.
# For now, avoiding object-oriented stuff.

# A model is represented by a list with at least:
#   vars - a named list of variables (the names should
#     be unique)
#   factors - a list of factors
# Each factor has a list consisting of:
#   m - the messages from the factor
#   update - function which takes a named list of messages,
#     and returns a named list of messages.
#   (later, maybe add log-evidence.)



# Adds in a message. The variable names in y should
# be a subset of the names in x.
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

# Constructs a factor with some set of methods.
# Args:
#   x - the variables connected to this factor
#   f - the "update" message for this factor
make.factor = function(x, f) {
  # the messages start at the belief, but zeroed out
  z = lapply(x, function(y) 0 * y )
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
  
  m = list(a2, a1)
  names(m) = names(a)
  m
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

