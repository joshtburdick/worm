# Computes coverage from BAM files.

# data.dir = "src/seq/coverage/lanes1and4/"

compute.coverage = function(data.dir, output.file) {
  x = NULL
  for(f in list.files(data.dir)) {
    cat(f, "")
    counts = read.table(gzfile(paste(data.dir, f, sep="")), sep="\t", as.is=TRUE)
    x = cbind(x, counts[,7])
    rownames(x) = counts[,4]
    colnames(x)[ncol(x)] = sub(".coverage", "", f)
  }
  write.table(x, file=gzfile(output.file),
    sep="\t", row.names=TRUE, col.names=NA)
}

# compute.coverage("git/unmix/seq/quant/Murray_050912/", "git/unmix/seq/quant/rawCoverage.tsv.gz")
# compute.coverage("src/seq/coverage/Murray_050912_as/", "R/unmix/sort_paper/seq/quant/rawCoverage_as.tsv")

compute.coverage("git/unmix/seq/quant/Murray_52831_092812/", "git/unmix/seq/quant/rawCoverage_Murray_52831_092812.tsv.gz")



