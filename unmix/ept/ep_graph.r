# Attempt at a more general tool for EP.

# A model is represented by a list with at least:
#   var - a named list of variables (the names should
#     be unique)
#   factor - a list of factors

# A variable is a list with at least these elements:
#   b - the total messages to/from this variable (as a matrix)
#   log.evidence - function which computes log-evidence

# Makes a variable.
make.var = function(b, log.evidence) list(b = b)

# A factor is a list with:
#   x - list of names of vars to which it's connected
#   to.f, from.f - lists of messages to and from this factor
#   update - function which, given messages to this factor ("to.f"),
#     sends messages from this factor ("from.f")
#   log.evidence - computes log-evidence for this factor (not yet implemented)

# Constructs a factor with some set of methods.
make.factor = function(x, f) {
  # the messages start at the belief, zeroed out
  z = lapply(x, function(a) 0 * a$b )

  list(x = x, to.f = z, from.f = z,
    update = f$update, log.evidence = f$log.evidence)
}

# An equality factor (maybe.)
eq.f = function(log.partition) function(a) {
  a1 = a[[1]]
  a2 = a[[2]]
  list(update = list(a2, a1),
    log.evidence = log.partition(a1+a2) - (log.partition(a1) + log.partition(a2)))
}


# Factor which constrains a variable (almost the same thing.)
constrain.f = function(log.partition, x) function(a) {
  a1 = a[[1]]

  list(update = list(x),
    log.evidence = log.partition(x+a1) - (log.partition(x) + log.partition(a1)))
}





# FIXME not working
# Updates a model some number of times, using "parallel" updating.
# Currently, no convergence testing.
#   m - a model, as a list with elements "vars" and "factors"
#   num.iters - number of iterations to do
# Returns: m (after some number of updates)
update.ep = function(m, num.iters = 30) {
  log.evidence = NULL

  for(iter in 1:num.iters) {

    # get messages from each variable
    sapply(m$factors, function(f)
      for(i in 1:length(f$x)) {
# cat("factor ", i, "\n")
# cat("dim(to.f) =", dim(f$to.f[[i]]), "\n")
# cat("dim(f$x...) =", dim(f$x[[i]]$b), "\n")
# cat("dim(from.f) =", dim(f$from.f[[i]]), "\n")
        f$to.f[[i]] = f$x[[i]]$b - f$from.f[[i]]
      })

    # clear the variables
    sapply(m$vars, function(x) if (!x$observed) x$b = x$b - x$b)

    # update messages from each factor
    sapply(m$factors, function(f) {
      print(f$from.f)
      print(f$to.f)
      print(f$update)
      f$from.f = f$update(f$to.f)
    })
cat("updated\n")
    # send messages to each variable
    sapply(m$factors, function(f)
      for(i in 1:length(f$x))
        if (!f$x[[i]]$observed)
          f$x[[i]]$b = f$x[[i]]$b + f$to[[i]])

    # compute log-evidence
    # FIXME: need this for variables as well
    le = sum(sapply(m$factors,
      function(f) f$log.evidence(f$to.f)), na.rm=TRUE)
    log.evidence = c(log.evidence, le)
  }

  m
}

