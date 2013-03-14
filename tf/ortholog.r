# Gets orthologs for TFs from biomaRt.
# XXX not working.

library("biomaRt")

# Biomart connections
ens.ce = useMart("ensembl", dataset="celegans_gene_ensembl")
ens.hs = useMart("ensembl", dataset="hsapiens_gene_ensembl")
ens.mm = useMart("ensembl", dataset="mmusculus_gene_ensembl")

# Generalized "conversion" function.
# Args:
#   mart - a biomaRt object
#   from - attribute to filter on
#   to - attribute to return
# Returns: function which converts names.
mart.convert = function(mart, from, to) function(x) {
  r = getBM(attributes = c(to),
    filters = c(from),   
    values = list(x),
    mart = mart)
  unique(r)
}

# Useful conversion functions.
wb.gene.name.to.gene.id =
  mart.convert(wb.gene, "gene_name", "sequence_name")
gene.id.to.entrez = 
  mart.convert(ensembl, "ensembl_gene_id", "entrezgene")



ens.attributes = listAttributes(ensembl)[,1]

a = ens.attributes[c(423:434)]

foo = getBM(attributes = c("ensembl_gene_id",
    "hsapiens_homolog_ensembl_gene",
    "hsapiens_homolog_perc_id_r1",
    "mmusculus_homolog_ensembl_gene",
    "mmusculus_homolog_perc_id_r1"),
  filters = c("wormbase_locus"),
  values = c("F38A6.1"),
  mart = ensembl.ce)


ens = useMart("ensembl", dataset="hsapiens_gene_ensembl")
if (FALSE) {
wb.gene = useMart("WS220", "wormbase_gene")

# this provides a whole bunch of matches
foo = getBM(
  attributes = c("public_name", "homolog_protein", "homolog_species", "homolog_info_accession"),
  filters = c("public_name"),
  values = c("pha-4"),
  mart = wb.gene)
}





