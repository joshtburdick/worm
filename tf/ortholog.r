# Gets orthologs for TFs from biomaRt.

library("biomaRt")

# list of TFs from wTF 2.1
wtf = read.table("data/tf/wTF2.1.csv", sep=",", header=TRUE, as.is=TRUE)
colnames(wtf) = c("sequence.name", "gene", "DNA.binding.domain",
  "array.coord", "pooling.coord", "X", "X.1", "clone.available",
  "source", "why.not.in.ORFeome", "fw.primer", "rev.primer",
  "gene.model", "gene.model.attempted.to.clone")
rownames(wtf) = wtf$sequence.name

# Biomart connections
ens.ce = useMart("ensembl", dataset="celegans_gene_ensembl")
ens.hs = useMart("ensembl", dataset="hsapiens_gene_ensembl")
ens.mm = useMart("ensembl", dataset="mmusculus_gene_ensembl")
ens.dm = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")

wb.gene = useMart("WS220", "wormbase_gene")

# MEME TF names
meme.tf = read.table("git/tf/motif/meme.tf.annotate.tsv",
  sep="\t", header=TRUE, as.is=TRUE)

# Generalized "conversion" function.
# Args:
#   mart - a biomaRt object
#   from - attribute to filter on
#   to - attribute to return
# Returns: function which converts names.
mart.convert = function(mart, from, to) function(x) {
  r = getBM(attributes = c(from, to),
    filters = c(from),   
    values = list(x),
    mart = mart)
  r
  r1 = r[ match(x, r[,1]) , 2 ]
  r1
}

# Useful conversion functions.
wb.gene.name.to.gene.id =
  mart.convert(wb.gene, "gene_name", "sequence_name")
gene.id.to.entrez =
  mart.convert(ensembl, "ensembl_gene_id", "entrezgene")

# Gets info. about whether there are known motifs.
# Args:
#   r - a data frame with columns "homolog_species"
#     and "homolog_gene"
# Returns: that data frame, with added column:
#   has_motif - 1 iff there is a correspondingly-named motif
#   motif - one of the motif names, chosen arbitrarily
#     (FIXME possibly this should include all of them)
check.for.motif = function(r) {
  r$has_motif = ""
  r$motif = NA

  r[ r$homolog_species == "Homo sapiens" &
    r$homolog_gene %in% meme.tf[ meme.tf$organism=="Hs" , "gene" ], "has_motif"] = 1
  r[ r$homolog_species == "Mus musculus" &
    r$homolog_gene %in% meme.tf[ meme.tf$organism=="Mm" , "gene" ], "has_motif"] = 1
  r[ r$homolog_species == "Drosophila melanogaster" &
    r$homolog_gene %in% meme.tf[ meme.tf$organism=="Dm" , "gene" ], "has_motif"] = 1

  # tack on id of one of the motifs
  m = meme.tf[ meme.tf$organism == "Hs" , ]
  r[ r$homolog_species == "Homo sapiens", "motif" ] =
    m[ match(r[ r$homolog_species == "Homo sapiens", "homolog_gene" ], m$gene) , "id" ]
  m = meme.tf[ meme.tf$organism == "Mm" , ]
  r[ r$homolog_species == "Mus musculus", "motif" ] =
    m[ match(r[ r$homolog_species == "Mus musculus", "homolog_gene" ], m$gene) , "id" ]
  m = meme.tf[ meme.tf$organism == "Dm" , ]
  r[ r$homolog_species == "Drosophila melanogaster", "motif" ] =
    m[ match(r[ r$homolog_species == "Drosophila melanogaster", "homolog_gene" ], m$gene) , "id" ]

  r[ is.na(r$motif), "motif" ] = ""

  r
}

# One way to compute orthologs, using Ensembl.
write.tf.ortho = function() {

tf1 = getBM(attributes = c("ensembl_gene_id",
    "hsapiens_homolog_ensembl_gene",
    "hsapiens_homolog_perc_id_r1",
    "mmusculus_homolog_ensembl_gene",
    "mmusculus_homolog_perc_id_r1",
    "dmelanogaster_homolog_ensembl_gene",
    "dmelanogaster_homolog_perc_id_r1"),
  filters = c("ensembl_gene_id"),
  values = c(as.character(wtf[,"sequence.name"])),
  mart = ens.ce)

tf.ortho = data.frame(gene.id = tf1$ensembl_gene_id,
  gene = wtf[tf1$ensembl_gene_id,"gene"],
  tf.type = wtf[tf1$ensembl_gene_id,"DNA.binding.domain"],
  hsapiens_homolog =
    mart.convert(ens.hs, "ensembl_gene_id", "hgnc_symbol")
      (tf1$hsapiens_homolog_ensembl_gene),
  hsapiens_homolog_perc_id_r1 = tf1$hsapiens_homolog_perc_id_r1,
  mmusculus_homolog =
    mart.convert(ens.mm, "ensembl_gene_id", "mgi_symbol")
      (tf1$mmusculus_homolog_ensembl_gene),
  mmusculus_homolog_perc_id_r1 = tf1$mmusculus_homolog_perc_id_r1,
  dmelanogaster_homolog =
    mart.convert(ens.dm, "ensembl_gene_id", "flybasename_gene")
      (tf1$dmelanogaster_homolog_ensembl_gene),
  dmelanogaster_homolog_perc_id_r1 = tf1$dmelanogaster_homolog_perc_id_r1,
  stringsAsFactors=FALSE)

tf.ortho[is.na(tf.ortho)] = ""

# look up whether there are motifs with that gene name
tf.ortho$has.hs.meme.motif = ifelse(tf.ortho$hsapiens_homolog %in%
  meme.tf[ meme.tf$organism=="Hs", "gene" ], 1, 0)
tf.ortho$has.mm.meme.motif = ifelse(tf.ortho$mmusculus_homolog %in%
  meme.tf[ meme.tf$organism=="Mm", "gene" ], 1, 0)
tf.ortho$has.dm.meme.motif = ifelse(tf.ortho$dmelanogaster_homolog %in%
  meme.tf[ meme.tf$organism=="Dm", "gene" ], 1, 0)
tf.ortho$has.ce.meme.motif = ifelse(toupper(tf.ortho$gene) %in%
  meme.tf[ meme.tf$organism=="Ce", "gene" ], 1, 0)

write.table(tf.ortho, file="git/tf/ortholog.tsv",
  sep="\t", row.names=FALSE, col.names=TRUE)
}

# Writes orthologs based on Wormbase data.
write.tf.ortho.wb = function() {

wb.tf = getBM(
  attributes = c("public_name", "sequence_name",             
"homolog_species",                      
"homolog_info_database",                 
"homolog_info_accession",                   
"homolog_blastp_rank",                
"homolog_blastp_evalue"),
  filters = c("sequence_name"),
  values = wtf$sequence.name,
  mart = wb.gene)
wb.tf = wb.tf[ !is.na( wb.tf$homolog_blastp_evalue ) , ]
wb.tf = wb.tf[ order(wb.tf$public_name, wb.tf$homolog_blastp_evalue) , ]
wb.tf$homolog_gene = NA

wb.tf$type = wtf[wb.tf$sequence_name,"DNA.binding.domain"]

# tack on symbols
i = wb.tf$homolog_species == "Homo sapiens" &
  wb.tf$homolog_info_database == "ENSEMBL"
wb.tf[ i , "homolog_gene" ] =
  mart.convert(ens.hs, "ensembl_peptide_id", "hgnc_symbol")(
    wb.tf[ i , "homolog_info_accession" ] )

i = wb.tf$homolog_species == "Mus musculus" &
  wb.tf$homolog_info_database == "SW"
wb.tf[ i , "homolog_gene" ] =
  mart.convert(ens.mm, "uniprot_swissprot_accession", "mgi_symbol")(
    wb.tf[ i , "homolog_info_accession" ] )

i = wb.tf$homolog_species == "Drosophila melanogaster" &
  wb.tf$homolog_info_database == "FLYBASE"
wb.tf[ i , "homolog_gene" ] =
  mart.convert(ens.dm, "flybasecgid_gene", "flybasename_gene")(
    wb.tf[ i , "homolog_info_accession" ] )


wb.tf = wb.tf[ !is.na(wb.tf$homolog_gene) , ]

wb.tf = check.for.motif(wb.tf)

write.table(wb.tf, file="git/tf/wb.ortholog.tsv",
  sep="\t", row.names=FALSE, col.names=TRUE)

}


# write.tf.ortho()
write.tf.ortho.wb()


