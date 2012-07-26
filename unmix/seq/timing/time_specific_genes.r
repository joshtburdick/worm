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
w = log.ets / apply(log.ets, 2, sum)

t.mean = apply( t(w) * time.points, 2, sum)
t.s2 = apply( t(w) * (time.points^2), 2, sum)
t.sd = sqrt( t.s2 - (t.mean^2) )

# plot these


#sorted.by.time = log.ets[ order(cor(t(log.ets), 1:23)) , ]

#write.table(round(sorted.by.time, 1), file="git/unmix/seq/timing/expr.by.time.tsv", sep="\t",
#  row.names=TRUE, col.names=NA)




