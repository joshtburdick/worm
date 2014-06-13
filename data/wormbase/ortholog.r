# Orthology information from Wormbase.

source("git/utils.r")
source("git/data/biomart_utils.r")

# Biomart connections
ens.hs = useMart("ensembl", dataset="hsapiens_gene_ensembl")
ens.mm = useMart("ensembl", dataset="mmusculus_gene_ensembl")
ens.dm = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")

# updating this from WS240     
r = read.table(gzfile("data/wormbase/c_elegans.PRJNA13758.WS240.orthologs.txt.gz"),
  sep="\t", fill=TRUE, as.is=TRUE)

# drop last line, as it's an "empty record"
r = r[ 1:(nrow(r)-1) , ]

# XXX somewhat funky hack
# find starts of records
record.start = cummax( ifelse(r[,1] == "=", 1:nrow(r), 0) )

# tack on relevant gene name for each record
r = cbind(worm.gene = as.character(r[record.start+1, 2]), r)
# r$worm.gene = as.character(r[record.start+1, 2])
r$worm.gene = as.vector(r$worm.gene)

# remove record separators, and first line of each record
r = r[ - c(record.start, record.start+1) , ]

colnames(r) = c("worm.gene", "species", "ensembl.gene.id", "method")

r$ensembl.gene.id = sub("^ENSEMBL:", "", r$ensembl.gene.id)

# write.tsv(r, "git/data/wormbase/ortholog.tsv")

r$gene = NA

# convert IDs for some of these
i = r$species == "Homo sapiens"
r[ i, "gene" ] =
  mart.convert(ens.hs, "ensembl_peptide_id", "hgnc_symbol")(r[i, "ensembl.gene.id"])
i = r$species == "Mus musculus"
r[ i, "gene" ] =
  mart.convert(ens.mm, "ensembl_peptide_id", "mgi_symbol")(r[i, "ensembl.gene.id"])
i = r$species == "Drosophila melanogaster"
r[ i, "gene" ] =
  mart.convert(ens.dm, "ensembl_peptide_id", "flybasename_gene")(r[i, "ensembl.gene.id"])

# for now, only keeping cases in which we know the gene name
# (and thus, only three species)
r = r[ !is.na(r$gene) , ]

write.tsv(r, "git/data/wormbase/ortholog.tsv")

