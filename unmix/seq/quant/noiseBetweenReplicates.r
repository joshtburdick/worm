# Attempting to model sampling error.

# get the raw read counts
r = as.matrix(read.table(gzfile("git/unmix/seq/quant/rawCoverage.tsv.gz"),
  sep="\t", as.is=TRUE, header=TRUE, row.names=1))
r = r[ rownames(r) != "ribosomal_RNA" , ]

# Gets the total counts for a gene across several technical replicates.
total.counts = function(pattern) {
  r1 = r[ , grep(pattern, colnames(r)) ]
  apply(r1, 1, sum)
}

# Plots mean vs. variance for several samples for which we have replicates.
plot.mean.and.var = function(r, main) {
  # normalize to each group of experiments having 1e6 reads
  r1 = t( t(r) / (apply(r, 2, sum) / 1e6) )
  a = data.frame(mean = log2(1+apply(r1, 1, mean)),
    var = log2(1+apply(r1, 1, var)))
  slope = lm(a$var ~ 0 + a$mean)$coef[1]

  smoothScatter(a$mean, a$var, pch=20, cex=0.2,
    main=main,
    xlab=expression(log[2]("1 + mean(reads)")),
    ylab=expression(log[2]("1 + var(reads)")))

  abline(0,1, col="#808080", lwd=2)
  mtext(paste("slope =", round(slope, 3)), cex=0.8)
  abline(0, slope, col="darkblue", lwd=2)
}

plot.mean.and.var.all = function() {
  pdf("git/unmix/seq/quant/noiseBetweenReplicates.pdf", width=7.5, height=7.5)
  par(mfrow=c(2,2))
  plot.mean.and.var(cbind(total.counts("cndm1p8_19"), total.counts("cndm1_12_14"),
    total.counts("cndm1p1_4")), "cnd-1 (n=3)")
  plot.mean.and.var(cbind(total.counts("pham4p9_1"), total.counts("pham4p12_9")),
    "pha-4 (n=2)")
  plot.mean.and.var(cbind(total.counts("cndm1_ungated"), total.counts("pham4_ungated")),
    "unsorted (n=2)")
  plot.mean.and.var(cbind(total.counts("cndm1_singlets"), total.counts("pham4_singlets")),
    "singlets (n=2)")
  dev.off()
}

plot.mean.and.var.all()


