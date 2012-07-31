
# Adds an observation to a normal variable.
# Args:
#   m - mean of observation
#   sd - standard deviation (if this is a "soft" observation)
norm.observe = function(m, sd = 0) list(

  update = function(   ) {


  },

  log.evidence = function(   ) {


  }

}


# Adds an observation to a multivariate normal.
mvnorm.constrain = function(A, b, b.var) list(

  update = function(a) {




  },

  # this is basically the likelihood for this observation
  log.evidence = function(a) {
    0    # FIXME

  }
)



