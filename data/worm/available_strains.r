# (Very) approximate catalog of available worm strains.

source("git/utils.r")

# list of transgenes, from WormBase
wb.transgene = read.table(gzfile(
  "data/wormbase/transgene_WS220.txt.gz"),
  sep="\t", header=TRUE, as.is=TRUE)
# for now, we only keep things of type "Drives_Transgene",
# as they seem more likely to indicate a reporter of
# expression (as opposed to that protein being expressed)
wb.transgene = wb.transgene[
  wb.transgene$Transgene.Type == "Drives_Transgene" , ]

cat("number of transgenes =",
  length(unique(wb.transgene$Gene.Public.Name)), "\n")



# Get names of strains available from the CGC.
cgc.strain.list = read.table(pipe(paste(
  "gunzip -c data/worm/cgc_strains_celelist_2.txt.gz |",
  " fgrep \"Strain:\" | cut -f2 -d ':' | tr -d ' ' | sort | uniq"), "r"),
  as.is=TRUE)[,1]






