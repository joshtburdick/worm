# Writes out a table for submission to the SRA.

source("git/utils.r")

exp.name = read.tsv("git/unmix/seq/quant/experimentNames1.tsv")
exp.name = exp.name[ ! grepl("HS |RNAi", exp.name$name) , ]

# Gets all the experiment names in a run.
exp.names.in.run = function(r) {
  a = read.tsv(paste0(
    "git/unmix/seq/quant/rawCoverage/WS220_20140111/", r, ".tsv.gz"))
  unique(sub("^0[0-9]_", "", colnames(a)))
}

# label approximate date of each experiment
exp.name$coll.date = NA
exp.name[ rownames(exp.name) %in%
  exp.names.in.run("20110922"), "coll.date"] =
  "20110922"
exp.name[ rownames(exp.name) %in%
  exp.names.in.run("Murray050912"), "coll.date"] =
  "20120509"
exp.name[ rownames(exp.name) %in%
  exp.names.in.run("Murray092812"), "coll.date"] =
  "20120928"

# create identifiers which avoid some potentially
# problematic characters (such as "(+-)")
id = exp.name$name
id = gsub("\\.", "", id)
id = gsub("F21D59", "F21D5.9", id)   # XXX
id = gsub("\\(\\+\\)", "plus", id)
id = gsub("\\(\\-\\)", "minus", id)
id = gsub(" ", "_", id)

# save mapping to identifiers
names(id) = exp.name$name
write.tsv(id, "git/unmix/seq/quant/SRA_BioSample_ID.tsv")

# create the table
r = data.frame(
  "sample_name" = id,
  "sample_title" = exp.name$name,
  organism = "C. elegans",
  strain = "",
  "dev_stage" = "mixed-stage embryo",
  "sex" = "hermaphrodite",
  "tissue" = paste(exp.name$name, "FACS-sorted cells"),
  "collected_by" = "Travis Walton",
  "collection_date" = exp.name$coll.date,
  "genotype" = "",
  stringsAsFactors=FALSE)
r = r[ order(r$"collection_date", r$"sample_title") , ]

write.table(r, file="git/unmix/seq/quant/SRA_BioSample.tsv",
  sep="\t", na="", row.names=FALSE, col.names=TRUE, quote=FALSE)
# strain and genotype need to be filled in by hand

