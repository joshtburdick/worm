# Gets orthology information about TFs from WormBase.

source("git/data/biomart_utils.r")

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

# Gets orthologs based on Wormbase data.
get.tf.ortho.wb = function() {

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

  wb.tf
}

# XXX freezing this for now
# write.table(get.tf.ortho.wb(), file="git/data/homology/wormbase.tf.ortholog.tsv",
#   sep="\t", row.names=FALSE, col.names=TRUE)

