# Attempt at a more general library for EP.
# For now, avoiding object-oriented stuff.

# Constructs a model.
# Args:
#   vars - a named list of variables. These should be
#     in terms of "natural parameters", so that adding
#     them multiplies the actual density functions.
#     (You may want to initialize these somehow.)
#   factors - a named list of lists, specifying factors.
#     Each element of this should have a factor, followed
#     by the names of the variables on which it depends.
# Returns: an EP model (represented as a list.)
#   (The messages are initially flat; I don't know if
#   there's a reason to alter them, but you can.)
ep.model = function(vars, factors) {
  
  # messages (initially flat)
  m = NULL
  for(i in names(factors)) {
    f1 = factors[[i]]
    m1 = NULL
    for(j in 2:length(f1)) {
      m1[[ f1 [[j]] ]] = 0 * vars[[ f1[[j]] ]]
    }
    m[[i]] = m1
  }

  # factors
  f = NULL
  for(i in names(f)) {
    f[[i]] = factors[[i]][[1]]
  }

  list(vars = vars, messages = m, factors = f)
}

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

# Updates a model once, using "parallel" updating.
# Currently, no convergence testing.
#   m - a model, as constructed above
# Returns: m, after one round of parallel updates
update.ep = function(m) {

  m1 = list()

  # loop through names of factors
  for(i in names(m$messages)) {

    # message to that factor
    m.to = message.subtract(m$vars, m$messages[[i]])

    # compute factor value
    m.f = m$factors[[i]](m.to)

    # compute message from this factor
    m1[[i]] = message.subtract(m.f, m$messages[[i]])
  }

  # clear the variables
  v1 = lapply(m$vars, function(v) 0 * v)

  # sum messages from each factor
  for(i in 1:length(f1)) {
    v1 = message.add(v1, m1[[i]])
  }

  list(vars = v1, messages = m1, factors = m$factors)
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
constrain.f = function(a) function(x) {
  x[[1]] = a
  x
}

