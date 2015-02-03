# Attempt to annotate the motifs in MEME with
# "organism" and "gene symbol" in a vaguely-consistent way.

r = read.table("data/tf/meme/motifList.tsv",
  sep="\t", header=TRUE, as.is=TRUE)

r$organism = ""
r$gene = ""

# a crude test for whether a gene symbol looks like a mouse gene
is.mouse.gene = function(x) grepl("^[A-Z][a-z]", x)

# first, include any motifs that had exactly the same
# name as an ortholog.
# XXX using previous version of this
wb.ortholog = read.table(
  "/home/jburdick/gcb/git/data/homology/wormbase.tf.ortholog.tsv",
  sep="\t", header=TRUE, as.is=TRUE)

i = r$name %in% wb.ortholog[
    wb.ortholog$homolog_species=="Homo sapiens", "homolog_gene" ]
r[ i, "gene" ] = r[ i, "name" ]
r[ i, "organism" ] = "Hs"

i = r$name %in% wb.ortholog[
    wb.ortholog$homolog_species=="Mus musculus", "homolog_gene" ]
r[ i, "gene" ] = r[ i, "name" ]
r[ i, "organism" ] = "Mm"

i = r$name %in% wb.ortholog[
    wb.ortholog$homolog_species=="Drosophila melanogaster", "homolog_gene" ]
r[ i, "gene" ] = r[ i, "name" ]
r[ i, "organism" ] = "Dm"


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



# getting more annotation from JASPAR database, somewhat hackily
matrix = read.table("data/tf/jaspar/MATRIX.txt",
    sep="\t", header=FALSE, as.is=TRUE)
colnames(matrix) =  c("motif.id", "motif.collection",
  "id1", "version", "gene")
matrix$id = paste0(matrix$id1, ".", matrix$version)
rownames(matrix) = matrix$id
matrix = matrix[ !duplicated( matrix$id ) , ]

# omitting species info in some cases for now
if (FALSE) {
  matrix.species = read.table("data/tf/jaspar/MATRIX_SPECIES.txt",
      sep="\t", header=FALSE, as.is=TRUE)
  colnames(matrix.species) = c("motif.id", "species.id")
  matrix1 = merge(matrix, matrix.species)
  matrix1$species = ""
  matrix1[ matrix1$species.id=="7227", "species" ] = "Dm"
  # matrix1[ matrix1$species.id=="10116", "species" ] = "Hs"
  # matrix1[ matrix1$species.id=="10090", "species" ] = "Mm"
}

i = r$id %in% rownames(matrix)
r[ i, "gene" ] = matrix[ r[i,"id"] , "gene" ]

r = r[ r$gene != "" , ]
r = r[ !duplicated(r$id) , ]

# write.table(r, file="git/tf/motif/meme.tf.annotate.tsv",
#   sep="\t", row.names=FALSE, col.names=TRUE)

