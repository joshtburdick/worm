# Runs all the tests and simulations.

# simulated datasets
source("sim.data/sim.cell.cor.r")
source("sim.data/sim.expr.cell.r")
source("sim.data/sim.expr.sym.r")

# pick reporters
source("pick_reporters.r")

# unmix using various methods
source("pseudoinverse/unmix.r")
source("EP/unmix.r")
source("trunc.pseudoinverse/unmix.r")
source("trunc.pseudoinverse.cor/unmix.r")
source("sampling/unmix.r")

# small test of EP (which doesn't require the actual data)
source("EP/approx_region_test.r")


