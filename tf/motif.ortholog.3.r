# Like motif.ortholog.2, but merging in the Hughes 2014 data.

source("git/utils.r")
source("git/data/name_convert.r")

# information about orthologs
motif.ortholog = read.table("git/tf/motif.ortholog.2.tsv",
  as.is=TRUE, header=TRUE)
motif.filter = read.tsv("git/tf/motif/motifFilter.tsv")
# add representative motif name, from the motif clustering
motif.ortholog$canonical.motif =
  motif.filter[ motif.ortholog$motif, "canonical.name" ]
motif.ortholog = motif.ortholog[ !is.na(motif.ortholog$canonical.motif) , ]
# for each motif, list of potential orthologs
orthologs.by.motif = by(motif.ortholog$gene, motif.ortholog$canonical.motif,
  function(x) {
    x = as.character(x)
    nhrs = grep("nhr", x)
    if (length(x) >= 5 && length(nhrs) >= 5) {
      return(unique(c(grep("nhr", x, invert=TRUE, value=TRUE), paste(length(nhrs), "NHRs"))))
    }

    return(unique(as.character(x))) 
  }
)

m1 = motif.ortholog[ , c("gene", "canonical.motif", "motif", "method") ]
colnames(m1) = c("gene", "motif.id", "motif.name", "method")

# hughes.motif = read.table("data/tf/hughes/TF_Information_all_motifs.tsv",
#   sep="\t", header=TRUE, as.is=TRUE)

# later, may use newer annotation...
hughes.motif = read.table(gzfile("data/tf/hughes/Caenorhabditis_elegans_2014_12_02_11_30_am/TF_Information_all_motifs_plus.txt.gz"),
  sep="\t", header=TRUE, as.is=TRUE)

m2 = hughes.motif[ hughes.motif$Motif_ID != "." ,
  c("TF_Name", "Motif_ID", "DBID.1", "MSource_Identifier") ]
m2[,4] = "Hughes"
colnames(m2) = colnames(m1)

motif.ortholog.3 = unique(rbind(m1, m2))

# possibly slightly improved motif names
motif.ortholog.3$motif.name.1 = motif.ortholog.3$motif.name

meme.tf.annotate = read.table("git/tf/motif/meme.tf.annotate.tsv", sep="\t", header=TRUE, as.is=TRUE)
meme.tf.annotate =
  meme.tf.annotate[ !duplicated(meme.tf.annotate$id) , ]
meme.tf.annotate$name.1 =
  paste(meme.tf.annotate$organism, meme.tf.annotate$gene)
rownames(meme.tf.annotate) = meme.tf.annotate$id

i = motif.ortholog.3$motif.id %in% rownames(meme.tf.annotate)
motif.ortholog.3[ i , "motif.name.1" ] =
  meme.tf.annotate[ motif.ortholog.3[ i , "motif.id" ], "name.1" ]

# also rename Hughes motifs, when there's a more readable name
i = motif.ortholog.3$motif.name.1 %in% rownames(meme.tf.annotate)
motif.ortholog.3[ i , "motif.name.1" ] =
  meme.tf.annotate[ motif.ortholog.3[ i , "motif.name.1" ], "name.1" ]

motif.ortholog.3 = motif.ortholog.3[ order(motif.ortholog.3$gene) , ]
write.tsv(motif.ortholog.3, "git/tf/motif.ortholog.3.tsv")

