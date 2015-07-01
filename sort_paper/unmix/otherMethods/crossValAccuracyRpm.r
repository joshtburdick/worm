# Tests predictions by cross-validation (taken from
# the pseudoinverse code.)

library(corpcor)

source("git/utils.r")
load("R/lineage/tree_utils.Rdata")

# the sort matrix, with positive and negative fractions
source("git/sort_paper/unmix/sortMatrix.r")
# ??? tweak the ceh-6 & hlh-16 so that the numbers
# are 0 or 1 ? (currently one cell is at 0.5)
m1 = {
  m.plus = m.unnormalized[1:14,]
  rownames(m.plus) = paste(rownames(m.plus), "(+)")
  # XXX omitting (-) sort fractions for hlh-16 & irx-1
  m.minus = 1 - m.unnormalized[c(1:5,8:14),]
  rownames(m.minus) = paste(rownames(m.minus), "(-)")

  # for now, omitting double-sorted fractions
  m1 = rbind(m.plus, m.minus)   # , m.unnormalized[15:18,])

  # move things away from 0/1?
#  m1 = 0.9 * m1 + 0.05

  m1
}
# m1 = m1 / apply(m1, 1, sum)

# the read data
rpm = as.matrix(read.tsv("git/cluster/readsPerMillion.tsv"))
rpm = rpm[ , !grepl("^(HS|RNAi)", colnames(rpm)) ]

# add on replicate averages (and remove replicates)
rpm = {
  avg = function(f) apply(rpm[ , f  ], 1, mean)

  cnd1.f = c("cnd-1 8/19", "cnd-1 12/14", "cnd-1 1/4")
  pha4.f = c("pha-4 5/9", "pha-4 9/1", "pha-4 12/9")
  replicates = paste(c(cnd1.f, pha4.f), c("(+)", "(-)"))

  # tack on averages
  rpm.avg = cbind(
    "cnd-1 (+)" = avg(paste(cnd1.f, "(+)")),
    "cnd-1 (-)" = avg(paste(cnd1.f, "(-)")),
    "pha-4 (+)" = avg(paste(pha4.f, "(+)")),
    "pha-4 (-)" = avg(paste(pha4.f, "(-)")))

  # remove replicates (and double negatives)
  rpm = rpm[ , setdiff(colnames(rpm), replicates) ]
  rpm = rpm[ , setdiff(colnames(rpm), "ceh-6 (-) hlh-16 (-)") ]
  cbind(rpm, rpm.avg)
}

# just using the genes with highest max expression
# in some sort fraction
max.expr = sort(apply(rpm, 1, max), decreasing=TRUE)

# for now, only using 5,000 most highly expressed;
# this leaves some for cross-validation
rpm = rpm[ names(max.expr)[ 1:15000 ] , ]

# tack on sum of other genes (in rpm)
x.other = 1e6 - apply(rpm, 2, sum)
rpm = rbind(rpm, "_other" = x.other)

# Definition of "enrichment."
log2.enrich = function(x.plus, x.minus)
  log2( 3 + x.plus ) - log2( 3 + x.minus )

# Plots accuracy, using cross-validation.
# Args:
#   unmix.f - an unmixing function (see write.crossval.graphs())
#   m - sort matrix
#   r - the read data, as reads per million
#   test.set - the fractions to use in testing
# Returns: list, indexed by test set, each containing a list of:
#   predicted.x.f - the predicted expression in that fraction
#   corr - correlation of predictions with the measured expression
#     of each gene
# Side effects: plots a scatterplot comparing these
plot.crossval.accuracy = function(unmix.f, m, r, test.set) {
  result = NULL

  mn = m / apply(m, 1, sum)

  for(a in test.set) {
    write.status(a)

    # test fractions to omit
    tf = paste(a, c("(+)", "(-)"))
    all.fractions = intersect(rownames(m), colnames(r))
    p = setdiff(all.fractions, tf)
    unmix.r = unmix.f(m[p,], r[,p])
    x = unmix.r$x

    measured.enrich = log2.enrich(r[, tf[1] ], r[ , tf[2] ])
    predicted.enrich =
      as.vector(log2.enrich(mn[ tf[1] , ] %*% t(x), mn[ tf[2] , ] %*% t(x)))
    corr = cor(measured.enrich, predicted.enrich)
    lim = range(c(measured.enrich, predicted.enrich))
    plot(measured.enrich, predicted.enrich,
      main=a, xlab="Measured enrichment", ylab="Predicted enrichment",
      xlim=lim, ylim=lim,
      pch=183, font=5, cex=0.5, col="#00000080", xaxt="n", yaxt="n")

    # XXX this is just to get the axis labels in the same font
    axis(1); axis(2) 
    mtext(paste("Correlation =", round(corr, 3)),
      cex=0.75, line=0.5)
    abline(0, 1, col="#00000080")
    result[[a]] = list(predicted.x = x, corr = corr)
  }

  result
}

# Computes cross-validation accuracy, and graphs comparisons,
# as a PNG image.
# Args:
#   unmix.f - an unmixing function, with args:
#     m - the sort matrix (with rows normalized)
#     r - the read ratios (as reads per million)
#   which returns a list including x, a matrix of predicted rpm.
#   output.name - what to call the output image
# Side effects: writes graphs, and correlation accuracies.
write.crossval.graphs = function(unmix.f, output.name) {
  output.dir = "git/sort_paper/unmix/otherMethods/"
  png(paste0(output.dir, "/", output.name, ".png"),
    width=1000, height=800)
  par(mfrow=c(4,5))

  a = plot.crossval.accuracy(unmix.f, m1, rpm,
    setdiff(single.fractions, c("hlh-16", "irx-1")))

  dev.off()

  correlations = sapply(a, function(x) x$corr)
  cat("correlations range = ", range(correlations), "\n")
  cat("median = ", median(correlations), "\n")
  cat("mean = ", mean(correlations), "\n")
}



