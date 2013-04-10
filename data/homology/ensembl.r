# Munging of Ensembl ortholog information
# (from Ensmart; I presume this is Compara)
# not sure this is better than what I had from Ensembl before

data.path = "data/homology/ensembl/"

source("git/data/biomart_utils.r")

# Biomart connections
ens.ce = useMart("ensembl", dataset="celegans_gene_ensembl")
ens.hs = useMart("ensembl", dataset="hsapiens_gene_ensembl")
ens.mm = useMart("ensembl", dataset="mmusculus_gene_ensembl")
ens.dm = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
wb.gene = useMart("WS220", "wormbase_gene")

read.homolog.file = function(f, organism) {
  r = read.table(gzfile(paste(
    data.path, f, sep="/")),
    sep="\t", header=TRUE, check.names=FALSE, as.is=TRUE)

  r1 = data.frame(gene=r[,1], ortholog=r[,3],
    type=r[,6],
    organism = organism,
    ancestor=r[,7],
    id.rel.to.gene = r[,8],
    id.rel.to.ortholog = r[,9], stringsAsFactors=FALSE)

  r1
}

r.hs = read.homolog.file("Ce_Hs_homologs_20130323.txt.gz", "H.sapiens")
r.mm = read.homolog.file("Ce_Mm_homologs_20130323.txt.gz", "M.musculus")
r.dm = read.homolog.file("Ce_Dm_homologs_20130323.txt.gz", "D.melanogaster")

# convert gene names
r.hs[,"ortholog"] =
  mart.convert(ens.hs, "ensembl_gene_id", "hgnc_symbol")(r.hs[,"ortholog"])
r.mm[,"ortholog"] =
  mart.convert(ens.mm, "ensembl_gene_id", "mgi_symbol")(r.mm[,"ortholog"])
r.dm[,"ortholog"] =
  mart.convert(ens.dm, "ensembl_gene_id", "flybasename_gene")(r.dm[,"ortholog"])

ens.ortholog = rbind(r.hs, r.mm, r.dm)
ens.ortholog$gene =
  mart.convert(wb.gene, "sequence_name", "public_name")(ens.ortholog$gene)

ens.ortholog = ens.ortholog[ !is.na(ens.ortholog$gene) , ]
ens.ortholog = ens.ortholog[ !is.na(ens.ortholog$ortholog) , ]


