# Attempt to annotate the motifs in MEME with
# "organism" and "gene symbol" in a vaguely-consistent way.

r = read.table("data/tf/meme/motifList.tsv",
  sep="\t", header=TRUE, as.is=TRUE)

r$organism = ""
r$gene = ""

# a crude test for whether a gene symbol looks like a mouse gene
is.mouse.gene = function(x) grepl("^[A-Z][a-z]", x)

r[ r$database %in% c("chen2008", "wei2010_human_mws") ,
  "organism" ] = "Hs"
i = r$database == "chen2008"
r[ i , "gene" ] = r[ i, "id" ]
i = r$database == "wei2010_human_mws"
r[ i , "gene" ] = sub("^h-", "", r[ i, "id" ])

r[ r$database %in% c("homeodomain", "uniprobe_mouse",
    "wei2010_mouse_mws", "wei2010_mouse_pbm") ,
  "organism" ] = "Mm"
i = r$database %in% c("homeodomain", "uniprobe_mouse")
r[ i, "gene" ] = sub("_.*", "", r[ i, "name" ])
i = r$database %in% c("wei2010_mouse_mws", "wei2010_mouse_pbm")
r[ i, "gene" ] = sub("^m-", "", r[ i, "id" ])

# omitting this, as genes from all organisms are jumbled together
# i = r$database == "zhao2011"
# r[ i, "gene" ] = sub("_.*", "", sub("^[^_]+_", "", r[ i, "id" ]))

r[ r$database %in% c("fly_factor_survey", "flyreg.v2") ,
  "organism" ] = "Dm"
i = r$database == "fly_factor_survey"
r[ i, "gene" ] = sub("_.*", "", r[ i, "name" ])
i = r$database == "flyreg.v2"
r[ i, "gene" ] = r[ i, "id" ]

r[ r$database %in% c("uniprobe_worm") ,
  "organism" ] = "Ce"
i = r$database == "uniprobe_worm"
r[ i, "gene" ] = r[ i, "name" ]

i = r$database == "jolma2013"
r[ i, "gene" ] = sub("_.*", "", r[ i, "id" ])
r[ i, "organism" ] = ifelse(is.mouse.gene(r[ i, "gene" ]), "Mm", "Hs")

r = r[ r$organism != "" & r$gene != "" , ]

write.table(r, file="git/tf/motif/meme.tf.annotate.tsv",
  sep="\t", row.names=FALSE, col.names=TRUE)



