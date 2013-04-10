# Summarizes ortholog information for Ce genes,
# compared to three other organisms, from
# four organisms.

source("git/utils.r")
source("git/data/biomart_utils.r")

# Biomart connections
ens.ce = useMart("ensembl", dataset="celegans_gene_ensembl")
ens.hs = useMart("ensembl", dataset="hsapiens_gene_ensembl")
ens.mm = useMart("ensembl", dataset="mmusculus_gene_ensembl")
ens.dm = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
wb.gene = useMart("WS220", "wormbase_gene")

# Gets Ensembl orthologs.
source("git/data/homology/ensembl.r")
r.e = data.frame(gene = ens.ortholog$gene,
  method = "Ensembl",
  species = ens.ortholog$organism,
  related.gene = ens.ortholog$ortholog,
  gene.score = ens.ortholog$id.rel.to.gene / 100,
  related.gene.score = ens.ortholog$id.rel.to.ortholog / 100,
  ortho.cluster = NA,
  type = ens.ortholog$type)

# Gets HomoloGene orthologs.
source("git/data/homology/homologene.r")

# Gets Homologene orthologs for one organism.
get.homologene.ortho = function(tax.id, species) {
  h.ce =
    homologene[ homologene[,"taxonomy.id"] == "Caenorhabditis elegans", ]
  h.other = homologene[ homologene[,"taxonomy.id"] == tax.id, ]
  i = match(h.ce$group.id, h.other$group.id)
  r.h = data.frame(gene = h.ce[,"gene"],
    method = "HomoloGene",
    species = species,
    related.gene = h.other[i,"gene"],
    gene.score = NA, related.gene.score = NA,
    ortho.cluster = h.ce$group.id, type = NA)
}
r.h = rbind(
  get.homologene.ortho("Homo sapiens", "H.sapiens"),
  get.homologene.ortho("Mus musculus", "M.musculus"),
  get.homologene.ortho("Drosophila melanogaster", "D.melanogaster"))
r.h = r.h[ !is.na(r.h$related.gene) , ]

# Gets InParanoid orthologs for one organism.
# First, read table, and rename genes.
in.paranoid = read.table(gzfile("git/data/homology/InParanoid.tsv.gz"),
  sep="\t", header=TRUE, as.is=TRUE)
in.paranoid[ in.paranoid$species=="C.elegans", "gene" ] =
  mart.convert(wb.gene, "gene", "public_name")(in.paranoid[in.paranoid$species=="C.elegans", "gene"])
in.paranoid[ in.paranoid$species=="H.sapiens", "gene" ] =
  mart.convert(ens.hs, "ensembl_gene_id", "hgnc_symbol")(in.paranoid[in.paranoid$species=="H.sapiens", "gene"])
in.paranoid[ in.paranoid$species=="M.musculus", "gene" ] =
  mart.convert(ens.mm, "ensembl_gene_id", "mgi_symbol")(in.paranoid[in.paranoid$species=="M.musculus", "gene"])
in.paranoid[ in.paranoid$species=="D.melanogaster", "gene" ] =
  mart.convert(ens.dm, "ensembl_gene_id", "flybasename_gene")(in.paranoid[in.paranoid$species=="D.melanogaster", "gene"])
in.paranoid = in.paranoid[ !is.na(in.paranoid$gene) , ]

# Reformats InParanoid orthologs for one organism.
get.in.paranoid.ortho = function(file, species) {
  in.paranoid.ce = in.paranoid[ in.paranoid$file==file &
    in.paranoid$species=="C.elegans", ]
  in.paranoid.other = in.paranoid[ in.paranoid$file==file &
    in.paranoid$species==species, ]
  i = match(in.paranoid.ce$cluster, in.paranoid.other$cluster)
  r.i = data.frame(gene=in.paranoid.ce[,"gene"],
    method="InParanoid",
    species = species,
    related.gene = in.paranoid.other[i,"gene"],
    gene.score = in.paranoid.ce[,"score"],
    related.gene.score = in.paranoid.other[i,"score"],
    ortho.cluster = in.paranoid.ce[,"cluster"],
    type=NA)
}
r.i = rbind(
  get.in.paranoid.ortho("InParanoid.C.elegans-H.sapiens.xml", "H.sapiens"),
  get.in.paranoid.ortho("InParanoid.C.elegans-M.musculus.xml", "M.musculus"),
  get.in.paranoid.ortho("InParanoid.C.elegans-D.melanogaster.xml", "D.melanogaster"))


# Get OrthoMCL orthologs.
ortho.mcl = read.table(gzfile("git/data/homology/ortho.mcl.tsv.gz"),
  sep="\t", row.names=1, header=TRUE, as.is=TRUE, quote="")

i = ortho.mcl[,"hs.gene"] != ""
r.o = data.frame(gene=ortho.mcl[i,"ce.gene"],
  method="OrthoMCL", species="H.sapiens",
  related.gene=ortho.mcl[i,"hs.gene"],
  gene.score=NA, related.gene.score=NA, ortho.cluster=NA, type=NA)

i = ortho.mcl[,"mm.gene"] != ""
r.o = rbind(r.o, data.frame(gene=ortho.mcl[i,"ce.gene"], 
  method="OrthoMCL",
  species="M.musculus",
  related.gene=ortho.mcl[i,"mm.gene"],
  gene.score=NA, related.gene.score=NA, ortho.cluster=NA, type=NA))

i = ortho.mcl[,"dm.gene"] != ""
r.o = rbind(r.o, data.frame(gene=ortho.mcl[i,"ce.gene"],
  method="OrthoMCL",
  species="D.melanogaster",
  related.gene=ortho.mcl[i,"dm.gene"],
  gene.score=NA, related.gene.score=NA, ortho.cluster=NA, type=NA))

# Gets data from previous WormBase query.
wb = read.table("git/data/homology/wormbase.tf.ortholog.tsv",
  sep="\t", header=TRUE)
r.w = data.frame(gene = wb$public_name, method = "WormBase",
  species = wb$homolog_species, related.gene = wb$homolog_gene,
  gene.score = NA, related.gene.score = -log10(wb$homolog_blastp_evalue),
  ortho.cluster = NA, type = wb$type)

ortho.five.methods = rbind(r.e, r.h, r.i, r.o, r.w)

ortho.five.methods = ortho.five.methods[
  !is.na(ortho.five.methods$gene) & !is.na(ortho.five.methods$related.gene) , ]
ortho.five.methods = unique(ortho.five.methods)

write.table(ortho.five.methods,
  file=gzfile("git/data/homology/ortho.five.methods.tsv.gz"),
  sep="\t", na="", row.names=FALSE, col.names=TRUE, quote=FALSE)

