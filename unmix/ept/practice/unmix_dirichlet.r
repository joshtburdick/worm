# Attempt at unmixing, modelling proportions as Dirichlet.

# Conditions a Dirichlet distribution, on some elements
# adding up to something.
# Args:
#   a - initial Dirichlet parameters
#   A. b - these give the linear constraint
# Returns: a counts, conditional on A(mean(Dirichlet(a))) = b
dirichlet.cond = function(a, A, b) {

  # the concentration
  s = apply(a, 1, sum)

  # the current sums
  m1 = apply(a*A, 1, sum)
  m2 = s - m1
  
  # what those sums should be
  s1 = s * b
  s2 = s * (1-b)

  # rescale (this should preserve the concentration)
  a1 = a*A*(s1/m1) + a*(1-A)*(s2/m2)

  a1[ is.na(a1) ] = 0
  a1
}

# Does unmixing.
#   A, b - these define an exact constraint.
#     A should be a 0/1 matrix (which is admittedly restrictive.)
unmix.dirichlet = function(A, b) {

  # the messages from each factor; initially flat
  m = 0 * A    # + 1

  # the posterior
  x = apply(m, 2, sum)

  for (iter in 1:100) {

    # compute messages to each factor
    m.to = t( x - t(m) )
    # XXX
#    for(i1 in 1:nrow(m.to))
#      m.to[i1,] = x

    # compute posterior "from" each factor
    mm = dirichlet.cond(m.to + 1, A, b) - 1
#    x1[x1<0] = 0

    # update messages from each factor
#    m.change = t( t(x1) - x )
    m = 0.5 * t( t(mm) - x ) + m

    # update posterior
    x = apply(m, 2, sum)

#    print(x)
  }

  list(m = m, x = x + 1, mean= (x+1)/sum(x+1))
}



A = rbind(
  c(1,1,1,1,1,0),
  c(0,0,0,1,1,1))

a1 = dirichlet.cond(a, A, c(0.5, 0.1))

r = unmix.dirichlet(A, c(0.5, 0.1))


