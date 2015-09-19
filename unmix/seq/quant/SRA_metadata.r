# Metadata for the read files in the SRA submission.

source("git/utils.r")

s = read.tsv("git/unmix/seq/quant/SRA_BioSample.tsv")

e = read.tsv("git/unmix/seq/quant/experimentNames1.tsv")
e = e[ ! grepl("^(HS |RNAi)", e$name) , ]

# create identifiers which avoid some potentially
# problematic characters (such as "(+-)")
id = e$name
id = gsub("\\.", "", id)
id = gsub("F21D59", "F21D5.9", id)   # XXX
id = gsub("\\(\\+\\)", "plus", id)
id = gsub("\\(\\-\\)", "minus", id)
id = gsub(" ", "_", id)
e$id = id

# the sample names, as submitted to the SRA
e$sample_name = sub("_rep_[0-9]", "", e$id)
e$sample_name[ e$sample_name=="lowInput" ] = "pha-4_ungated"
e$sample_name[ e$sample_name=="normalInput" ] = "pha-4_ungated"

# mapping from sample name to SRA sample ID
samn = {
  samn1 = read.tsv("git/unmix/seq/quant/BioSampleObjects.tsv")
  samn = rownames(samn1)
  names(samn) = samn1[,1]
  samn
}


# which experiments were newer (paired-end and amplified)
i = !(e$id %in% c("lowInput", "normalInput",
  "pha-4_plus_rep_1", "pha-4_minus_rep_1"))

r = data.frame(bioproject_accession = "PRJNA295677",
  sample_name = samn[ e$sample_name ],
  library_ID = e$id,
#  title = "RNA-Seq of FACS-sorted embryonic C. elegans cells",
  library_strategy = "RNA-Seq",
  library_source = "TRANSCRIPTOMIC",
  library_selection = ifelse(i, "other", "RANDOM"),
  library_layout = ifelse(i, "Paired", "Single"),

  platform = "ABI_SOLID",
  instrument_model = "AB 5500xl Genetic Analyzer",
  design_description = ifelse(i,
    "Embryonic cells were FACS-sorted by marker expression. RNA was extracted using a RNeasy kit (Qiagen), and poly-A RNA was amplified using a T7 RNA polymerase aRNA protocol (Ambion MessageAMP II aRNA kit).",
    "Embryonic cells were FACS-sorted by marker expression; RNA was extracted using a RNeasy kit (Qiagen)."),
  filetype = "BAM",
  # these will be filled in later
  filename = "", filename2 = "", filename3 = "",
  check.names=FALSE, stringsAsFactors=FALSE)

# Get mapping from original filenames to filenames in submission.
# This assumes that all the experiments for a given run are in a
# particular subdirectory (which currently is true.)
p = "/murrlab/seq/tophat2/WS220_20140111/"

# mapping from files to experiments
file.to.experiment = NULL

# loop through the directories of .bam files
for(p1 in c("20110922", "Murray050912", "Murray092812")) {
  files = list.files(paste0(p, p1))
  files = sort(grep(".bai", files, value=TRUE, invert=TRUE))
  file.sample = sub("^0[0-9]_", "", sub("\\.bam$", "", files))
  
  # scan through the experiments
  for(i in 1:nrow(e)) {
    e1 = rownames(e)[i]
    f1 = files[ file.sample == e1 ]
    if (length(f1) > 0) {
      new.filenames = paste0(paste(r[i,"library_ID"], 1:length(f1),
        sep="_"), ".bam")
      r[i, c(12:(11+length(f1)))] = new.filenames

      file.to.experiment = rbind(file.to.experiment,
        data.frame(f = paste0(p, p1, "/", f1),
          new.f = new.filenames,
          stringsAsFactors = FALSE))
    }
  }
}

write.tsv(file.to.experiment, "git/unmix/seq/quant/file.to.experiment.tsv")

write.table(r, "git/unmix/seq/quant/SRA_metadata.tsv",
  sep="\t", row.names=FALSE)



