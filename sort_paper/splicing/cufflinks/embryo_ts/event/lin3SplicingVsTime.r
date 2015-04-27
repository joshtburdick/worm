# Plots the aforementioned timeseries.

source("git/utils.r")

r = read.tsv("git/sort_paper/splicing/cufflinks/embryo_ts/event/cufflinksSummary.tsv")
colnames(r)[1] = "time"
r = r[ order(r$time) , ]

pdf("git/sort_paper/splicing/cufflinks/embryo_ts/event/lin3SplicingVsTime.pdf")

ylim = range(c(r$FPKM_conf_lo, r$FPKM, r$FPKM_conf_hi))

plot(0, 0, type="n",
  xlim=range(r$time), ylim=ylim,
  main="LIN-3 alternative splicing",
  xlab="Time (minutes)", ylab="Expression (FPKM)")

cols = hsv(c(0,1,2,3)/4, 1, 0.9, 0.8)

# XXX not actually genes
genes = sort(unique(r$Gene))

for(j in c(1:4)) {
  s = genes[j]
  par(new=TRUE)
  r1 = r[ r$Gene == s , ]
  plot(r1$time, r1$FPKM, ylim=ylim, type="l", col = cols[j],
    main="", xlab="", ylab="", xaxt="n", yaxt="n", lwd=2)
  arrows(r1$time, r1$FPKM_conf_lo, r1$time, r1$FPKM_conf_hi,
    length=0.07,
    col=cols[j], code=3, angle=90)
}

legend("topright", legend=genes, col=cols, lwd=2)

dev.off()

