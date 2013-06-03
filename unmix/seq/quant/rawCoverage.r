# Computes coverage from BAM files.

# data.dir = "src/seq/coverage/lanes1and4/"
base.path = "git/unmix/seq/quant/rawCoverage/"

compute.coverage = function(data.dir, output.file) {
  data.dir = paste(base.path, data.dir, sep="/")
  output.file = paste(base.path, output.file, sep="/")

  x = NULL
  for(f in list.files(data.dir)) {
    cat(f, "")
    counts = read.table(gzfile(paste(data.dir, f, sep="")), sep="\t", as.is=TRUE)
    x = cbind(x, counts[,7])
    rownames(x) = counts[,4]
    colnames(x)[ncol(x)] = sub(".tsv.gz$", "", f)
  }
  write.table(x, file=gzfile(output.file),
    sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
}

if (FALSE) {
compute.coverage("20110922/", "rawCoverage_20110922.tsv.gz")
compute.coverage("20110922_as/", "rawCoverage_20110922_as.tsv.gz")

compute.coverage("Murray_050912/", "rawCoverage_Murray_050912.tsv.gz")
compute.coverage("Murray_050912_as/", "rawCoverage_Murray_050912_as.tsv.gz")

compute.coverage("Murray_52831_092812/", "rawCoverage_Murray_52831_092812.tsv.gz")
compute.coverage("Murray_52831_092812_as/", "rawCoverage_Murray_52831_092812_as.tsv.gz")
}

compute.coverage("embryo_ts/", "rawCoverage_embryo_ts.tsv.gz")

