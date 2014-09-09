# Attempt at a more general library for EP.

# A model is represented by a list with at least:
#   vars - a named list of variables (the names should
#     be unique)
#   factors - a list of factors

# Makes a model.
# Args:
#   vars - list of variables
#   factors - list of factors
make.model = function(vars, factors) {
  # XXX this admittedly looks pretty silly
  r = list(vars = vars, factors = factors)
#  attr(r, "class") = "ep.model"
  r
}

# A variable is a list with at least these elements:
#   b - the total messages to this variable (as a matrix)
#   observed - TRUE iff this variable is observed
#   log.evidence - function which computes log-evidence

# Makes a variable.
make.var = function(b, log.evidence) {
  r = list(b = b, observed = FALSE, log.evidence = log.evidence)
#  attr(r, "class") = "ep.var"
  r
}

# A factor is a list with:
#   from.f - list of messages from this factor
#   update - function which, given messages to this factor,
#     sends messages from this factor
#   log.evidence - function which computes log-evidence for this factor

# Constructs a factor with some set of methods.
# Args:
#   x - the variables connected to this factor
#   f - list containing the "update" and "log.evidence" functions
#     for this factor
make.factor = function(x, f) {
  # the messages start at the belief, but zeroed out
  z = lapply(x, function(a) 0 * a$b )
  names(z) = names(x)

  list(from.f = z,
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
# ??? this may be all wrong.
constrain.f = function(log.partition, x) function(a) {
  a1 = a[[1]]

  list(update = list(x),
    log.evidence = log.partition(x+a1) - (log.partition(x) + log.partition(a1)))
}






# FIXME not working
# Updates a model once, using "parallel" updating.
# Currently, no convergence testing.
#   m - a model, as a list with elements "vars" and "factors"
# Returns: m, after one round of parallel updates
update.ep = function(m) {


  # get messages from each factor
  f1 = sapply(m$factors, function(f)
    for(i in 1:length(f$x)) {
# cat("factor ", i, "\n")
# cat("dim(to.f) =", dim(f$to.f[[i]]), "\n")
# cat("dim(f$x...) =", dim(f$x[[i]]$b), "\n")
# cat("dim(from.f) =", dim(f$from.f[[i]]), "\n")
        f$to.f[[i]] = f$x[[i]]$b - f$from.f[[i]]
      })

  # clear the variables
  v1 = lapply(m$vars, function(v) { v$b = 0 * v$b; v })

  # send messages from each factor



  list(vars = v1, factors = f1)
}

