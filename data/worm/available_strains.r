# (Very) approximate catalog of available worm strains.

source("git/utils.r")

wb.transgene = read.table(gzfile(
  "data/worm/Wormbase_strain_drives_transgene.tsv.gz"),
  sep="\t", as.is=TRUE)
colnames(wb.transgene) =
  c("name", "reporter", "gene", "strain", "method")
wb.transgene = wb.transgene[ wb.transgene$gene != "" , ]
wb.transgene = wb.transgene[ wb.transgene$strain != "" , ]

# Get names of strains available from the CGC.
# XXX hack
cgc.strain.list = read.table(pipe(paste(
  "gunzip -c data/worm/cgc_strains_celelist_2.txt.gz |",
  " fgrep \"Strain:\" | cut -f2 -d ':' | tr -d ' ' | sort | uniq"), "r"),
  as.is=TRUE)[,1]

wb.transgene$at.CGC = wb.transgene$strain %in% cgc.strain.list

wb.transgene = unique(wb.transgene[,c(2,3,4,6)])

# summarize this per-gene
wb.transgene$summary = paste0(wb.transgene$strain,
  "(", wb.transgene$reporter,
  ifelse(wb.transgene$at.CGC, ",at CGC", ""), ")")

wb.transgene = wb.transgene[ order(wb.transgene$at.CGC, decreasing=TRUE) , ]

wb.transgene.by.gene = c(by(wb.transgene$summary,
  wb.transgene$gene,
  function(x) paste(x, collapse=" ")))
wb.transgene.by.gene = sapply(wb.transgene.by.gene,
  function(x) if (nchar(x) >= 80) paste0(substr(x,1,80), "...") else x)




