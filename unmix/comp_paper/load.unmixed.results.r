# Loads assorted data sets for comparing unmixing accuracy,
# reading from .Rdata objects.

p = "~/gcb/git/unmix/unmix_comp/src/"

expr.cell = as.matrix(read.table(
  paste(p, "../data/exprCell.tsv", sep="/"),
  sep="\t", header=TRUE, row.names=1, as.is=TRUE))

load(paste(p, "sim.data/sim.cell.cor.Rdata", sep="/"))
load(paste(p, "sim.data/sim.expr.cell.Rdata", sep="/"))
load(paste(p, "sim.data/sim.expr.sym.Rdata", sep="/"))

source("R/unmix/comp_paper/expression.on.off.r")

num.reporters = c(10,20,30,50,75,100)

# Loads one set of data.
get.unmixed = function(base.name) {
  r = list()

  for (dataset in c("expr.cell", "sim.cell.cor", "sim.expr.cell", "sim.expr.sym")) {
    r[[dataset]] = list()

    for(nr in num.reporters) {
      unmix.result = NULL
      load(paste(p, base.name, "/", dataset, ".", nr, ".Rdata", sep=""))
      if (!is.null(unmix.result)) {
        r[[dataset]][[nr]] = unmix.result$x
      }
    }
  }

  r
}

# Loads in all the data.
unmix.r = list(
  pseudo = get.unmixed("pseudoinverse"),
  tp = get.unmixed("trunc.pseudoinverse"),
  tpc = get.unmixed("trunc.pseudoinverse.cor"),
  ep = get.unmixed("EP"),
  ep2 = get.unmixed("EP.2"),
  ln = get.unmixed("lognormal"),
  mf = get.unmixed("mult.fractions"),
  nnls = get.unmixed("nnls")
)

# definition of when genes are "on" or "off" (somewhat hokey)
expr.on.off = list(expr.cell = expr.cell >= 2000,
  sim.cell.cor = sim.cell.cor >= 0.2,
  sim.expr.cell = sim.expr.cell >= 5,
  sim.expr.sym = sim.expr.sym >= 5)

# NB: this is all deprecated
if (FALSE) {
# load the pseudoinverse
p.expr.cell = get.unmixed("pseudoinverse/expr.cell")
p.sim.cell.cor = get.unmixed("pseudoinverse/sim.cell.cor")
p.sim.expr.cell = get.unmixed("pseudoinverse/sim.expr.cell")
p.sim.expr.sym = get.unmixed("pseudoinverse/sim.expr.sym")

# load the EP predictions
ep.expr.cell = get.unmixed("EP/expr.cell")
ep.sim.cell.cor = get.unmixed("EP/sim.cell.cor")
ep.sim.expr.cell = get.unmixed("EP/sim.expr.cell")
ep.sim.expr.sym = get.unmixed("EP/sim.expr.sym")

# load the truncated pseudoinverse
tp.expr.cell = get.unmixed("trunc.pseudoinverse/expr.cell")
tp.sim.cell.cor = get.unmixed("pseudoinverse/sim.cell.cor")
tp.sim.expr.cell = get.unmixed("trunc.pseudoinverse/sim.expr.cell")
tp.sim.expr.sym = get.unmixed("trunc.pseudoinverse/sim.expr.sym")

# load the truncated pseudoinverse with correlation
tpc.expr.cell = get.unmixed("trunc.pseudoinverse.cor/expr.cell")
tpc.sim.cell.cor = get.unmixed("trunc.pseudoinverse.cor/sim.cell.cor")
tpc.sim.expr.cell = get.unmixed("trunc.pseudoinverse.cor/sim.expr.cell")
tpc.sim.expr.sym = get.unmixed("trunc.pseudoinverse.cor/sim.expr.sym")
}

