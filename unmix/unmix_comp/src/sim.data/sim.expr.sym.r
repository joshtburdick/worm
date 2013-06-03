# Creates a data set with symmetric lineages on.

load("../data/tree_utils.Rdata")
source("sim.data/sim.expr.r")
source("sim.data/sim.expr.gamma.r")

# definitions of which lineages are symmetric (from Sulston)
sym.lineages = read.csv("../data/symmetric.lineages.csv", as.is=TRUE)

# Given a lineage, converts it to a symmetric lineage (if possible.)
sym.lineage = function(x) {
#  cat(x, "")
  p = function(a) paste("^", a, sep="")

  for(i in 1:nrow(sym.lineages)) {
    if (length(grep(p(sym.lineages[i,"lineage1"]), x)) > 0)
      return(sub(p(sym.lineages[i,"lineage1"]), sym.lineages[i,"lineage2"], x))

    if (length(grep(p(sym.lineages[i,"lineage1"]), x)) > 0)
      return(sub(p(sym.lineages[i,"lineage1"]), sym.lineages[i,"lineage2"], x))
  }

  return(NA)
}

# compute list of lineages with at least some number of cells on
cell.lineage.1 = cell.lineage.matrix[apply(cell.lineage.matrix, 1, sum)>=5,]

# compute all "big enough" pairs of symmetric lineages
sym.pairs = data.frame(lin = rownames(cell.lineage.1), stringsAsFactors=FALSE)
sym.pairs$sym.lin = sapply(sym.pairs$lin, sym.lineage)
sym.pairs = sym.pairs[ !is.na(sym.pairs$sym.lin) , ]
sym.pairs = sym.pairs[ sym.pairs$lin %in% rownames(cell.lineage.1) , ]
sym.pairs = sym.pairs[ sym.pairs$sym.lin %in% rownames(cell.lineage.1) , ]

cell.pairs.on = cell.lineage.1[ sym.pairs$lin , ] + cell.lineage.1[ sym.pairs$sym.lin , ]

set.seed(42)
sim.expr.sym = sim.expr(cell.pairs.on)

save(sim.expr.sym, file="sim.data/sim.expr.sym.Rdata")

# gamma-distributed version of this
sim.expr.sym.gamma = sim.expr.gamma(cell.pairs.on)
save(sim.expr.sym.gamma, file="sim.data/sim.expr.sym.gamma.Rdata")

