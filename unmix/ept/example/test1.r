# Some example EP models.

source("git/unmix/ept/ep_graph.r")
source("git/unmix/ept/normal.r")

x1 = norm.var(3)



m1 = list(vars = list(x1),
  factors = list(pf = make.factor(list(a = x1), positive.factor)))
    


