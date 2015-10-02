# Messages for EP, defined using S4 classes.
# Draws on Hadley Wickham's S4 tutorial, at
# http://adv-r.had.co.nz/S4.html

# the Message class

# FIXME arguably this should check that x is a list
# of numeric things, or Messages (or possibly lists)
check_Message = function(object) TRUE
setClass("Message", representation(x = "list"),
  validity = check_Message)

# Lift some one-var functions.
# (Omitting Math for now, as they don't seem useful for EP.)
setMethod("Summary", signature(c(x="Message")),
  function(x) {
    r = c()
    for(a in names(x@x))
      r = c(r, callGeneric(x@x[[a]]))
    0
    callGeneric(r)
  }
)

# Lift binary operations.
# Binary operations, numeric %op% Message...
setMethod("Ops", signature(e1="numeric", e2="Message"),
  function(e1, e2) {
    for(a in names(e2@x)) {
      e2@x[[a]] = callGeneric(e1, e2@x[[a]])
    }
# XXX unfortunately, this doesn't work
#    e2@x = lapply(x, function(y) callGeneric(e1, y))
    validObject(e2)
    return(e2)
  }
)

# ...and Message %op% numeric.
setMethod("Ops", signature(e1="Message", e2="numeric"),
  function(e1, e2) {
    for(a in names(e1@x)) {
      e1@x[[a]] = callGeneric(e1@x[[a]], e2)
    }
    validObject(e1)
    return(e1)
  }
)

# Lastly, Message %op% Message.
setMethod("Ops", signature(e1="Message", e2="Message"),
  function(e1, e2) {
    r = list()

    # if both elements are present, apply the op
    for(a in intersect(names(e1@x), names(e2@x)))     
      r[[a]] = callGeneric(e1@x[[a]], e2@x[[a]])

    # any element that's missing is treated as 0
    for(a in setdiff(names(e1@x), names(e2@x)))     
      r[[a]] = callGeneric(e1@x[[a]], 0)
    for(a in setdiff(names(e2@x), names(e1@x)))     
      r[[a]] = callGeneric(0, e2@x[[a]])

    e1@x = r
    validObject(e1)
    return(e1)
  }
)

