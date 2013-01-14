# Definition of when a gene is "on".

# determine when expression is "on" or "off" (defined using logistic regression)
expression.on.off = function(x, int.off=1000, int.on=3000, int.sd=1000) {
  p.off = dnorm(x, mean=int.off, sd=int.sd, log=TRUE)
  p.on = dnorm(x, mean=int.on, sd=int.sd, log=TRUE)

  # normalize
  a = p.on - p.off       # log-odds-ratio that this is on
  is.on = exp(a) / (exp(a) + exp(-a))
  is.on[ a > 100 ] = 1
  is.on[ a < -100 ] = 0

  is.on
}


