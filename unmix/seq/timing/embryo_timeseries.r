# Summarizes the embryonic timeseries data.

data.path = "git/unmix/seq/timing/embryo_ts_coverage/"

r = NULL

for(f in list.files("git/unmix/seq/timing/embryo_ts_coverage")) {
  cat(f, "")
  sample.name = sub(".tsv.bz2", "", f)
  x = read.table(paste(data.path, "/", f, sep=""),
    sep="\t", header=TRUE, as.is=TRUE, row.names=1)
  r = cbind(r, x[,3])
  rownames(r) = rownames(x)
  colnames(r)[ncol(r)] = sample.name
}

embryo.timeseries = r
write.table(r, file=gzfile("git/unmix/seq/timing/embryo.timeseries.tsv.gz"),
  sep="\t", row.names=TRUE, col.names=NA)

# normalize to "ppm"
embryo.timeseries.rpm =
  t( t(embryo.timeseries) / (apply(embryo.timeseries,2,sum) / 1e6) )

# estimate mean and sd of when each gene is on
time.dist = {

#  x1 = embryo.timeseries.rpm[ apply(embryo.timeseries.rpm, 1, mean) >= 1 , ]
  w = embryo.timeseries.rpm / apply(embryo.timeseries.rpm, 1, sum)
  w[ is.na(w) ] = 0

  t.mean = apply( t(w) * time.points, 2, sum)
  t.s2 = apply( t(w) * (time.points^2), 2, sum)
  t.sd = sqrt( t.s2 - (t.mean^2) )
  cbind(mean = t.mean, sd = t.sd)
}

# Writes out all genes sorted by time.
sorted.by.time =
  embryo.timeseries.rpm[ order(time.dist[,"mean"]) , ]
write.table(round(sorted.by.time, 4),
  file="git/unmix/seq/timing/expr.rpm.sorted.by.time.tsv", sep="\t",
    row.names=TRUE, col.names=NA)



