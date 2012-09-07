# Plots scatterplots between replicates.

# get read counts 
r = as.matrix(read.table("git/unmix/seq/quant/readsPerMillion.tsv",
  sep="\t", as.is=TRUE, header=TRUE, row.names=1, check.names=FALSE))
r = r[ rownames(r) != "ribosomal_RNA" , ]

lr = log2(1+r)

# Plots a scatterplot of two samples.
plot.scatter = function(x, y, xlab, ylab) {
  lim = c(0, max(x, y))

#  xlab.1 = paste("1 + mean(", xlab, " reads)", sep="")
#  ylab.1 = paste("1 + mean(", ylab, " reads)", sep="")

  plot(x, y, cex=0.1, pch=20, xlim=lim, ylim=lim,
    xlab=paste(xlab, "log2 coverage"),
    ylab=paste(ylab, "log2 coverage"))
# FIXME these labels aren't working
#    xlab=expression(log[2](xlab.1)),
#    ylab=expression(log[2](ylab.1)))
}

plot.scatter.all = function() {
  pdf("git/unmix/seq/quant/replicateScatterPlot.pdf", width=7.5, height=7.5)
  par(mfrow=c(2,2))

  plot.scatter(lr[,"cnd-1p8_19"], lr[,"cnd-1_12_14"],
    "cnd-1 (8/19)", "cnd-1 (12/14)")
  plot.scatter(lr[,"pha-4p9_1"], lr[,"pha-4p12_9"],
    "pha-4 (9/1)", "pha-4 (12/9)")
  plot.scatter(lr[,"cnd-1_singlets"], lr[,"pha-4_singlets"],
    "cnd-1 singlets", "pha-4 singlets")
  plot.scatter(lr[,"cnd-1_ungated"], lr[,"pha-4_ungated"],
    "cnd-1 ungated", "pha-4 ungated")

  dev.off()
}

plot.scatter.all()

