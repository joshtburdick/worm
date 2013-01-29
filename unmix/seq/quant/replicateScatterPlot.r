# Plots scatterplots between replicates.
# FIXME rename this, as it's no longer just about replicates.

# get read counts 
r = as.matrix(read.table("git/unmix/seq/quant/readsPerMillion.tsv",
  sep="\t", as.is=TRUE, header=TRUE, row.names=1, check.names=FALSE))
r = r[ rownames(r) != "ribosomal_RNA" , ]

lr = log2(1+r)

# singlet average
singlet.mean =
  (lr[,"cnd-1_singlets"] + lr[,"pha-4_singlets"]) / 2

# Plots a scatterplot of two samples.
plot.scatter = function(x, y, xlab, ylab) {
  lim = c(0, max(x, y))
  par(mar=c(5,4,4,2)+0.5)
#  xlab.1 = paste("1 + mean(", xlab, " reads)", sep="")
#  ylab.1 = paste("1 + mean(", ylab, " reads)", sep="")

  plot(x, y, cex=0.8, pch=20, xlim=lim, ylim=lim,
    xlab=paste(xlab, "log2 coverage"),
    ylab=paste(ylab, "log2 coverage"), cex.lab=2, cex.axis=2)
  abline(0,1, lwd=2)
# FIXME these labels aren't working
#    xlab=expression(log[2](xlab.1)),
#    ylab=expression(log[2](ylab.1)))

#  mtext(paste("r =", round(cor(x, y, use="pairwise.complete.obs"), 3)), cex=2)
}

plot.replicate.scatter = function() {
#  pdf("git/unmix/seq/quant/replicateScatterPlot.pdf", width=7.5, height=7.5)
  png("git/unmix/seq/quant/replicateScatterPlot.png", width=1000, height=1000)
  par(mfrow=c(2,2))
  par(mar=c(5,4,4,2)+0.5)
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

plot.scatterplots = function() {
  output.dir = "git/unmix/seq/quant/scatterplot/"
  system(paste("mkdir -p", output.dir))

  plot.1 = function(y, ylab) {
cat(ylab, "")
    png(paste(output.dir, y, ".png", sep=""),
      width=1000, height=1000)
    plot.scatter(singlet.mean, lr[,y], "singlets", ylab)
    dev.off()
  }

  for(g in c("ceh-6", "ceh-26", "unc-130", "hlh-16")) {
    plot.1(g, g)
  }

  plot.1("cnd-1_12_14", "cnd-1 (12-14)")

  png(paste(output.dir, "cnd-1_ungated.png", sep=""),
    width=1000, height=1000)
  plot.scatter(lr[,"cnd-1_singlets"], lr[,"cnd-1_ungated"],
    "cnd-1 singlets", "cnd-1 ungated")
  dev.off()

  png(paste(output.dir, "pha-4_ungated.png", sep=""),
    width=1000, height=1000)
  plot.scatter(lr[,"pha-4_singlets"], lr[,"pha-4_ungated"],
    "pha-4 singlets", "pha-4 ungated")
  dev.off()
}

plot.replicate.scatter()

plot.scatterplots()

