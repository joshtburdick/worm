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

colors = rainbow(length(facs.genes), v = 0.8, alpha=0.7)
names(colors) = facs.genes

pdf("git/sort_paper/FACS/timing/sortMarkerTimeProfile.pdf",
  width=8, height=6)

xlim = range(as.numeric(colnames(r.ts)))
ylim = range(r.ts.facs)

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

dev.off()

