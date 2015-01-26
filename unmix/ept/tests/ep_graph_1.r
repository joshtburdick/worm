# Simple test of EP graph code.

source("git/unmix/ept/ep_graph.r")
source("git/unmix/ept/gamma.r")


# Tiny model, with just one variable with a prior.
p = gamma.mv2n(rbind(m=c(1,2,3), v=c(1,1,1)))
vars = list(x1 = gamma.mv2n(rbind(m=c(1,1,1), v=c(1,1,1))))
vars = list(x1 = rbind(m=c(3,5,9), v=c(1,1,2.3456)))

m0 = ep.model(vars,
  list(f1 = list(prior.f(p), "x1"),
    f2 = list(prior.f(p), "x1")))




# A factor which constrains its args to sum to 1.
gamma.sum1.factor = function(a) {
  x = a[[1]]
  x1 = gamma.conditional.approx.1.beta(x)
  r = list(x1)
  names(r) = names(a)
  r
}

vars = list(x1 = gamma.mv2n(rbind(m=c(1,2,3), v=c(1,1,1))))

m1 = ep.model(vars, list(f1=list(gamma.sum1.factor, "x1")))

m1a = ep.update(m1)


