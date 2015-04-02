# Plots the aforementioned timeseries.

source("git/utils.r")

r = read.tsv("git/sort_paper/splicing/cufflinks/embryo_ts/event/groupedIsoforms.tsv")

times = as.numeric(rownames(r))

pdf("git/sort_paper/splicing/cufflinks/embryo_ts/event/lin3SplicingVsTime.pdf")

ylim = range(r)

plot(0, 0, type="n",
  xlim=range(times), ylim=ylim,
  main="LIN-3 isoform expression",
  xlab="Time (minutes)", ylab="Expression (FPKM)")

cols = hsv(c(0,1,2,3)/4, 0.7, 0.9)

for(j in c(1:4)) {
  par(new=TRUE)
  plot(times, r[,j], ylim=ylim, type="l", col = cols[j],
    main="", xlab="", ylab="", xaxt="n", yaxt="n", lwd=2)
}

legend("topright", legend=colnames(r), col=cols, lwd=2)

dev.off()

