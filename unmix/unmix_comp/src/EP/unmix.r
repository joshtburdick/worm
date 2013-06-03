# Does unmixing using EP.

library(corpcor)

source("EP/approx_region.r")

source("unmix_test.r")

# Does unmixing using EP.
unmix.ep = function(m, x.f) {

  # for now, normalizing this
  scale = max(x.f)
  x.f = x.f / scale

#  hack which seems to improve convergence
  eps = 1e-3

  x.f = as.vector(x.f + m %*% (rep(eps, ncol(m))))

  ep = approx.region(m, x.f, 0 * x.f,
    converge.tolerance = 1e-6, prior.var=Inf, max.iters = 200)

  x = ep$m - rep(eps, ncol(m))

  list(x = x, t = ep$t, update.stats = ep$update.stats,
    scale = scale, eps = eps,
    reporters = rownames(m), x.f = x.f)
}

# run.unmix(unmix.ep, "EP/")

