# Translates names of genes from OrthoMCL.

source("git/utils.r")
source("git/data/biomart_utils.r")

r = read.tsv(gzfile("git/data/homology/orthoMCL_genes_raw.txt.gz"))
rownames(r) = sub(":", "", rownames(r))

# Biomart connections
ens.ce = useMart("ensembl", dataset="celegans_gene_ensembl")
ens.hs = useMart("ensembl", dataset="hsapiens_gene_ensembl")
ens.mm = useMart("ensembl", dataset="mmusculus_gene_ensembl")
ens.dm = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
wb.gene = useMart("WS220", "wormbase_gene")

# FIXME: write this out in "pairwise" format
# translate names of genes
ortho.mcl = data.frame(
  ce.gene =
    mart.convert(wb.gene, "gene", "public_name")(r$cele),
  hs.gene =
    mart.convert(ens.hs, "ensembl_peptide_id", "hgnc_symbol")
      (r$hsap),
  mm.gene =
    mart.convert(ens.mm, "ensembl_peptide_id", "mgi_symbol")
      (r$mmus),
  dm.gene =
    mart.convert(ens.dm, "ensembl_peptide_id", "flybasename_gene")
      (r$dmel), stringsAsFactors=FALSE)

write.tsv(ortho.mcl, file=gzfile("git/data/homology/ortho.mcl.tsv.gz"))

