# Plots profile of sort markers w.r.t. time.

source("git/utils.r")

r = read.tsv("git/cluster/readRatios.tsv")
r.ts = as.matrix(r[ , c(24:38) ])
colnames(r.ts) =
  as.numeric(sub("^t.", "", colnames(r.ts)))

# XXX
facs.genes = c("ceh-26", "cnd-1", "pha-4", colnames(r)[7:18])
facs.genes = setdiff(facs.genes, c("mir-57", "pros-1"))
r.ts.facs = r.ts[ facs.genes , ]
rownames(r.ts.facs) = sub("ceh-26", "pros-1", rownames(r.ts.facs))

facs.genes = sort(sub("ceh-26", "pros-1", facs.genes))

# Yanai 2014 data
r.y = as.matrix(read.tsv(gzfile(
  "git/data/expr/hashimshony2014/whole.embryo.interval.timecourse.tsv.gz")))
r.y = r.y[ , c(-1,-2) ]
colnames(r.y) = sub("^mpfc_0?", "", colnames(r.y))
r.y = log2(1+r.y)

xlim=range(as.numeric(c(colnames(r.ts.facs), colnames(r.y))))
ylim.ts = range(r.ts.facs, na.rm=TRUE)
ylim.y = range(r.y)

# XXX for lin-3
r.ts.facs = r.ts

# Plots timeseries data from both datasets for one gene.
plot.ts = function(g) {

  plot(as.numeric(colnames(r.ts.facs)), r.ts.facs[g,],
    xlim=xlim, ylim=ylim.ts,
    main=g, xlab="time (minutes)", ylab="modENCODE timeseries expression",
    type="l", lwd=2, col="#000000a0")

  par(new=TRUE)
  plot(colnames(r.y), r.y[g,],
    xlim=xlim, ylim=ylim.y,
    main="", xlab="", ylab="", xaxt="n", yaxt="n",
    type="l", lwd=2, col="#0000ffa0")
  axis(4, col.axis="#0000ffff", col.ticks="#0000ffff")
  mtext("Yanai expression, log2(1+RPM)",
    side=4, line=3, col="#0000ffff", cex=0.75)
}

pdf("git/sort_paper/FACS/timing/lin-3.pdf",
  width=6, height=4)
#  width=7.5, height=9)
# par(mfrow=c(3,2))
par(mar=c(5,4,4,5)+0.1)

for(g in c("lin-3"))   # facs.genes)
  plot.ts(g)

# the original graph, with all the genes from the modENCODE timeseries
if (FALSE) {
  xlim = range(as.numeric(colnames(r.ts)))
  ylim = range(r.ts.facs)
  colors = rainbow(length(facs.genes), v = 0.8, alpha=0.7)
  names(colors) = facs.genes
  plot(0, 0, xlim=xlim, ylim=ylim, type="n",
    xlab="time (minutes)", ylab="relative expression")

  for(g in rownames(r.ts.facs)) {
    par(new=TRUE)
    plot(colnames(r.ts), r.ts.facs[ g , ],
      xlim=xlim, ylim=ylim, type="l",
      xlab="", ylab="", xaxt="n", yaxt="n",
      col=colors[g], lwd=3)
  }

  legend(250, -1, facs.genes, col=colors, cex=0.75, lwd=3)
}

dev.off()

