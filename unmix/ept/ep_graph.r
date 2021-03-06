# Attempt at a more general tool for EP.

# A variable is a list with at least two elements:
#   b - the total messages to/from this variable (as a matrix)
#   observed - boolean, indicating whether this is observed


# A factor is a list with:
#   x - list of variables to which it's connected
#   to.f, from.f - lists of messages to and from this factor
#   update - function which, given messages to this factor ("to.f"),
#     sends messages from this factor ("from.f")
#   log.evidence - computes log-evidence for this factor (not yet implemented)
# ??? should this include both messages "from"?

# Constructs a basic factor (this is a sort of "superclass constructor.")
new.factor = function(x) {
  z = sapply(x, function(a) a$b)
  z = z - z   # set these all to 0

  list(x = x, to.f = z, from.f = z)
}


# Updates a model some number of times. Currently, no convergence testing.
#   m - a model, as a list with elements "vars" and "factors"
#   num.iters - number of iterations to do
# Returns: m (after some number of updates)
update.ep = function(m, num.iters = 30) {


  for(iter in 1:num.iters) {

    # get messages from each variable
#    sapply(m$factors, functionf)

    # clear the variables
    sapply(m$vars, function(x) x$b = x$b - x$b)

    # send messages to each variable



    # FIXME: compute log-evidence

  }

  m
}

