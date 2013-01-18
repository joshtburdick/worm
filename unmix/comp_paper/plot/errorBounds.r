# Comparisons of error bounds.

library("ggplot2")

expr.cell = as.matrix(read.table("git/unmix/unmix_comp/data/exprCell.tsv",
  sep="\t", header=TRUE, row.names=1, as.is=TRUE))

load("git/unmix/comp_paper/mirror.sampling.summary.Rdata")
load("git/unmix/comp_paper/rda.sampling.summary.Rdata")
load("git/unmix/comp_paper/sampling/plot/samplingStats.Rdata")

load("git/unmix/unmix_comp/src/EP/lineage.totals.Rdata")


# Computes z-scores for sampling.
sampling.z = function(sampling.summary) {
  sampling.mu = sampling.summary[,"mean",]
  sampling.sd = sampling.summary[,"sd",]

  z1 = (expr.cell - sampling.mu) / sampling.sd
  z1 = z1[ !is.na(sampling.mu) ]

  # ??? is this valid?
  z1[ z1 <= -30 ] = -30
  z1[ z1 > 30 ] = 30
  z1
}


# Plots the histograms of z-scores.
plot.sampling.z.histograms = function() {
  pdf("git/unmix/comp_paper/plot/samplingZ.pdf",
    width=6, height=6)
#  par(mfrow=c(1,1))
  hist(sampling.z(samplingStats), breaks=100, col="grey",
    xlab="z",
    main="Expression, relative to sampling prediction")  
#  hist(sampling.z(mirror.sampling.summary), breaks=100, col="grey",
#    xlab="z-score of expression, relative to prediction",
#    main="Mirror sampling")  
  dev.off()
}

# Compares EP and sampling mean and variance.
ep.sampling.mean.var.plot = function() {

  png("git/unmix/comp_paper/plot/EP.and.sampling.png",
    width=1600, height=800)
  par(mfrow=c(1,2))
  par(mar=c(5.5,4.5,4,2) + 0.1)
  plot(samplingStats[,"mean",],
    lineage.totals[,,"per.cell.mean"],
    xlim=c(0,250000), ylim=c(0,250000),
    pch=20, col="#00000030", cex=2,
    main="Prediction mean",
    xlab="Random directions sampling", ylab="EP",
    cex.axis = 1.8, cex.lab=1.8, cex.main=2.2)
  abline(0, 1, col="#00000080", lwd=2)
  plot(samplingStats[,"sd",],
    sqrt(lineage.totals[,,"per.cell.var"]),
    xlim=c(0,250000), ylim=c(0,250000),
    pch=20, col="#00000030", cex=2,
    main="Prediction standard deviation",
    xlab="Random directions sampling", ylab="EP",
    cex.axis = 1.8, cex.lab=1.8, cex.main=2.2)
  abline(0, 1, col="#00000080", lwd=2)

  dev.off()
}

plot.sampling.z.histograms()
ep.sampling.mean.var.plot()

