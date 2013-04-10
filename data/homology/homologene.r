# Gets Homologene homology info.

r = read.table("data/homology/homologene.build67.data",
  sep="\t", fill=TRUE, as.is=TRUE, quote="")

colnames(r) = c("group.id", "taxonomy.id", "gene.id",
  "gene", "protein.gi", "protein")

# r = r[ nchar(r[,4]) <= 100 , ]

tax.ids.1 = read.table("data/homology/homologene.taxid_taxname.tsv",
  sep="\t", fill=TRUE, as.is=TRUE)
tax.ids = tax.ids.1[,2]
names(tax.ids) = tax.ids.1[,1]

r$taxonomy.id = tax.ids[ as.character( r$taxonomy.id ) ]

homologene = r



