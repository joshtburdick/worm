# Some example EP models.

source("git/unmix/ept/ep_graph.r")
source("git/unmix/ept/normal.r")
source("git/unmix/ept/positive_factor.r")

x1 = norm.var(3)

p = cbind(m = rep(0,3), v=rep(0,3))

m1 = list(vars = list(x1),
  factors = list(
    prior = make.factor(list(x = x1), normal.prior.factor(p)),
    pf = make.factor(list(x = x1), positive.factor)))
    


