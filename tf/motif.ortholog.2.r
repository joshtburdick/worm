# Gets orthologs for TFs from biomaRt.

library("biomaRt")

# ortholog annotation
ortho = read.table(gzfile(
# XXX not rebuilding this
  "/home/jburdick/gcb/git/data/homology/ortho.five.methods.tsv.gz"),
  sep="\t", header=TRUE, quote="")


# list of TFs from wTF 2.
wtf = read.table("data/tf/wTF2.1.csv", sep=",", header=TRUE, as.is=TRUE)
colnames(wtf) = c("sequence.name", "gene", "DNA.binding.domain",
  "array.coord", "pooling.coord", "X", "X.1", "clone.available",
  "source", "why.not.in.ORFeome", "fw.primer", "rev.primer",
  "gene.model", "gene.model.attempted.to.clone")
rownames(wtf) = wtf$sequence.name

# Biomart connections
if (TRUE) {
ens.ce = useMart("ensembl", dataset="celegans_gene_ensembl")
ens.hs = useMart("ensembl", dataset="hsapiens_gene_ensembl")
ens.mm = useMart("ensembl", dataset="mmusculus_gene_ensembl")
ens.dm = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")

wb.gene = useMart("WS220", "wormbase_gene")
}

# MEME TF names
meme.tf = read.table("git/tf/motif/meme.tf.annotate.tsv",
  sep="\t", header=TRUE, as.is=TRUE)

meme.tf.1 = data.frame(species = "",
  related.gene = meme.tf$gene,
  motif.id = meme.tf$id, stringsAsFactors=FALSE)
meme.tf.1[ meme.tf$organism=="Hs", "species" ] = "H.sapiens"
meme.tf.1[ meme.tf$organism=="Mm", "species" ] = "M.musculus"
meme.tf.1[ meme.tf$organism=="Dm", "species" ] = "D.melanogaster"

# Creates table of Ce genes and motif names. This should include both
# components of dimer motifs.
get.ce.motif.names = function() {

  # get all annotation for Ce genes
  r = meme.tf[ meme.tf$organism=="Ce", ]

  # separate into monomers and dimers
  dimers = c(grep("_", r$name, value=TRUE), "ttx-3::ceh-10", "MXL-1_MDL-1")
  r1 = r[ !(r$name %in% dimers) , ]
  r2a = r[ (r$name %in% dimers) , ]
  r2b = r2a
  # XXX add each half of dimers
  r2a[ , "gene" ] = "HLH-2"
  r2a[ r2a$name == "ttx-3::ceh-10" , "gene" ] = "ttx-3"
  r2a[ r2a$name == "MXL-1_MDL-1", "gene" ] = "mxl-1"
  r2b[ , "gene" ] = sub("HLH-2_", "", r2b$gene)
  r2b[ r2b$name == "ttx-3::ceh-10" , "gene" ] = "ceh-10"
  r2b[ r2b$name == "MXL-1_MDL-1", "gene" ] = "mdl-1"

  r3 = rbind(r1, r2a, r2b)
  z = data.frame(gene = tolower(r3$gene), method = "Ce_motif",
    species = "C.elegans", related.gene = tolower(r3$name),
    gene.score = NA, related.gene.score = NA,
    ortho.cluster = NA, type = "",
    has.motif = 1, motif = r3$id)
}

# Gets info. about whether there are known motifs.
# Args:
#   r - a data frame with columns "species"
#     and "related.gene"
# Returns: that data frame, with added column:
#   has_motif - 1 iff there is a correspondingly-named motif
#   motif - one of the motif names, chosen arbitrarily
#     (FIXME possibly this should include all of them)
check.for.motif = function(r) {
  r$has.motif = ""
  r$motif = NA

  r[ r$species == "H.sapiens" &
    r$related.gene %in% meme.tf[ meme.tf$organism=="Hs" , "gene" ], "has.motif"] = 1
  r[ r$species == "M.musculus" &
    r$related.gene %in% meme.tf[ meme.tf$organism=="Mm" , "gene" ], "has.motif"] = 1
  r[ r$species == "D.melanogaster" &
    r$related.gene %in% meme.tf[ meme.tf$organism=="Dm" , "gene" ], "has.motif"] = 1

  # tack on id of one of the motifs
  m = meme.tf[ meme.tf$organism == "Hs" , ]
  r[ r$species == "H.sapiens", "motif" ] =
    m[ match(r[ r$species == "H.sapiens", "related.gene" ], m$gene) , "id" ]
  m = meme.tf[ meme.tf$organism == "Mm" , ]
  r[ r$species == "M.musculus", "motif" ] =
    m[ match(r[ r$species == "M.musculus", "related.gene" ], m$gene) , "id" ]
  m = meme.tf[ meme.tf$organism == "Dm" , ]
  r[ r$species == "D.melanogaster", "motif" ] =
    m[ match(r[ r$species == "D.melanogaster", "related.gene" ], m$gene) , "id" ]

  r[ is.na(r$motif), "motif" ] = ""

  r
}

motif.ortholog = check.for.motif(ortho)
# add in Ce motifs
motif.ortholog = rbind(check.for.motif(ortho), get.ce.motif.names())
motif.ortholog.2 = motif.ortholog[ motif.ortholog$has.motif != "" , ]
motif.ortholog.2 =
  motif.ortholog.2[ order(motif.ortholog.2$gene, motif.ortholog.2$related.gene) , ]

write.table(motif.ortholog.2, file="git/tf/motif.ortholog.2.tsv",
  sep="\t", row.names=FALSE, col.names=TRUE)

