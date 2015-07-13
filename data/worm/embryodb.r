# Utilities for reading embryoDB.

library("XML")

# Given an EmbryoDB series name, gets the base name of the relevant directory.
embryodb.to.annot.dir = function(series) {
  x = xmlTreeParse(paste("/gpfs/fs0/l/murr/embryoDB/", series, ".xml", sep=""))
  
  as.vector(xmlRoot(x)["annots"][[1]]$attr)
}

# Given an EmbryoDB series name, gets one of the CSV files.
read.embryodb.dat.file = function(series, file.pattern="SCD") {
  annot.dir = embryodb.to.annot.dir(series)
# cat("annot.dir =", annot.dir, "\n")

  file.name = list.files(paste(annot.dir, "/dats", sep=""), pattern=paste("^", file.pattern, sep=""))
  f = paste(annot.dir, "dats", file.name, sep="/")

# cat("reading ", f, "\n")

  if (file.exists(f))
    read.csv(f, as.is=TRUE, header=TRUE)
  else
    NULL
}






