# Simple test of EP graph code.

source("git/unmix/ept/ep_graph.r")
source("git/unmix/ept/dirichlet.r")

a1 = c(1:5)
class(a1) = "dirichlet"
a2 = c(3,2,3,2,3)
class(a2) = "dirichlet"


vars = list(x1 = 0*a1, x2 = 0*a1, x3 = 0*a1)




