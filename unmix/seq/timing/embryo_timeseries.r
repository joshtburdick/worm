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



