# Plots accuracy of pseudoinverse using cross-validation.

source("git/sort_paper/unmix/pseudoinverseEnrichment.r")

output.dir = "git/sort_paper/unmix/pseudoinverseEnrichment/"

# Plots accuracy, using cross-validation.
# Args:
#   m - sort matrix
#   r - the read ratios
#   predict.set - the fractions to use in predicting
#   test.set - the fractions to use in testing
# Returns: list, indexed by test set, each containing a list of:
#   predicted.x.f - the predicted expression in that fraction
#   corr - correlation of predictions with the measured expression
# Side effects: plots a scatterplot comparing these
plot.crossval.accuracy = function(m, r, predict.set, test.set) {
  result = NULL

  for(a in test.set) {
    write.status(a)
    p = setdiff(predict.set, a)
# browser()
    x = t( pseudoinverse(m[p,]) %*% t(r[,p]) )
    predicted.x.f = as.vector( m[a,] %*% t(x) )
    corr = cor(r[,a], predicted.x.f)
    lim = range(c(r[,a], predicted.x.f))

    par(mar=c(4,4,4,2)+0.1)
    plot(r[,a], predicted.x.f,
      main=a, xlab="Measured enrichment", ylab="Predicted enrichment",
      xlim=lim, ylim=lim, cex.main=1.8, cex.lab=1.5, cex.lab=1.5,
      pch=183, font=5, cex=0.5, col="#00000080", xaxt="n", yaxt="n")

    # XXX this is just to get the axis labels in the same font
    axis(1, cex.axis=1.5); axis(2, cex.axis=1.5) 
    mtext(paste("Correlation =", round(corr, 3)),
      cex=0.85, line=0.1)
    abline(0, 1, col="#00000080")
    result[[a]] = list(predicted.x.f = predicted.x.f, corr = corr)
  }

  result
}

# Computes cross-validation accuracy, and graphs comparisons.
write.crossval.graphs = function() {
  pdf(paste0(output.dir, "/crossvalidationAccuracy.pdf"))

  a = plot.crossval.accuracy(m, r, single.fractions,
    c(single.fractions, "ceh-6 (+) hlh-16 (+)",
      "ceh-6 (+) hlh-16 (-)", "ceh-6 (-) hlh-16 (+)"))
  dev.off()

  correlations = sapply(a, function(x) x$corr)

  cat("correlations range = ", range(correlations), "\n")
  cat("median = ", median(correlations), "\n")
  cat("mean = ", mean(correlations), "\n")
}

# Same, but plots on one figure, as a png.
write.crossval.graphs.png = function() {
  png(paste0(output.dir, "/crossvalidationAccuracy.png"),
    width=800, height=1050)
#    width=1000, height=800)   was 4x5, with these dimensions
#  pdf(paste0(output.dir, "/crossvalidationAccuracy.pdf"),
#    width=14, height=12)
  par(mfrow=c(5,4))

  a = plot.crossval.accuracy(m, r, single.fractions,
    c(single.fractions, "ceh-6 (+) hlh-16 (+)",
      "ceh-6 (+) hlh-16 (-)", "ceh-6 (-) hlh-16 (+)")[1:5])

  dev.off()

  correlations = sapply(a, function(x) x$corr)

  cat("correlations range = ", range(correlations), "\n")
  cat("median = ", median(correlations), "\n")
  cat("mean = ", mean(correlations), "\n")
}

# write.crossval.graphs()
write.crossval.graphs.png()

