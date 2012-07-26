# Looks for genes which are only on at particular times.

embryo.timeseries = read.table("git/unmix/seq/timing/embryo.timeseries.tsv.gz",
  sep="\t", header=TRUE, row.names=1, as.is=TRUE)
time.points = c(0,30,60,90,120,140,180,240,270,300,330,360,390,
  420,450,480,540,570,600,630,660,690,720)
colnames(embryo.timeseries) = time.points

# normalize to "ppm"
embryo.timeseries =
  t( t(embryo.timeseries) / (apply(embryo.timeseries,2,sum) / 1e6) )

# log-transformed expression
log.ets = log2(1+embryo.timeseries)

# scatterplots one pair of genes.
plot.pair = function(i, j) {
  smoothScatter(log2(1+embryo.timeseries[,as.character(i)]),
    log2(1+embryo.timeseries[,as.character(j)]),
    xlab=paste("log2(1+coverage at time", i, ")"),
    ylab=paste("log2(1+coverage at time", j, ")"))
  abline(0,1, col="darkgrey")
}

# plots distribution of log-ratio of each pair of genes.
plot.all.log.ratios = function() {
  pdf("git/unmix/seq/timing/logRatioPlots.pdf", width=7.5, height=10)
#  par(mfrow=c(5,1))
  par(mfrow=c(3,2))
  for(i in time.points)
    for(j in time.points[ time.points > i ]) {
      cat(i, j, "   ")
      plot.pair(i, j)

      lr = log2(embryo.timeseries[,as.character(i)])
        - log2(embryo.timeseries[,as.character(j)])
      lr = lr[ !is.na(lr) ]
      lr = lr[ lr > -10 & lr < 10 ]
#      hist(lr, main = paste("log2 (", i, "/", j, ")", sep=""),
#        breaks=100, col="grey")
  }
  dev.off()
}

pdf("git/unmix/seq/timing/logRatioPlots1.pdf", width=9, height=3)
par(mfrow=c(1,3))
plot.pair(60,120)
plot.pair(390,630)
plot.pair(690,720)
dev.off()

# estimate mean and sd of when each gene is on
x1 = embryo.timeseries[ apply(embryo.timeseries, 1, mean) >= 10 , ]
w = x1 / apply(x1, 1, sum)
w[ is.na(w) ] = 0
time.dist = {
  t.mean = apply( t(w) * time.points, 2, sum)
  t.s2 = apply( t(w) * (time.points^2), 2, sum)
  t.sd = sqrt( t.s2 - (t.mean^2) )
  cbind(mean = t.mean, sd = t.sd)
}


# plot these
pdf("git/unmix/seq/timing/timeMeanAndSD.pdf", width=6, height=6)
smoothScatter(time.dist[,"mean"], time.dist[,"sd"],
  main="Time when genes are expressed",
  xlab="mean of time expressed", ylab="s.d. of time expressed")
print(quantile(time.dist[,"sd"], 0.2))
abline(h=quantile(time.dist[,"sd"], 0.2), col="#303030")
dev.off()

# Writes out all genes sorted by time. FIXME: use the above t.s2 and t.sd?
write.sorted.by.time = function() {
  #sorted.by.time = log.ets[ order(cor(t(log.ets), 1:23)) , ]
  #write.table(round(sorted.by.time, 1), file="git/unmix/seq/timing/expr.by.time.tsv", sep="\t",
  #  row.names=TRUE, col.names=NA)
}

# Plots expression of those genes.
r = as.matrix(read.table(gzfile("git/unmix/seq/quant/readsPerMillion.tsv.gz"),
  sep="\t", check.names=FALSE, header=TRUE, row.names=1, as.is=TRUE))

td1 = time.dist[ time.dist[,"sd"] <= quantile(time.dist[,"sd"], 0.2) , ]
td1 = td1[ order(td1[,"mean"], decreasing=FALSE) , ]

plot.time.profile = function(g, main) {
  plot(td1[,"mean"], log2(1 + r[ rownames(td1), g ]),
    main = main, xlab = "time", ylab="log2(1 + coverage)",
    pch=20, cex=0.5)
}
pdf("git/unmix/seq/timing/perStageMarkersInSorted.pdf", width=9, height=8)
par(mfrow=c(3,1))
plot.time.profile("all", "singlets")
plot.time.profile("pha-4", "pha-4")
plot.time.profile("ceh-26", "ceh-26")
dev.off()


