# Writes out a table of which genes have orthologs in several organisms.

source("git/utils.r")

ce.metaphors.gene = read.tsv("git/data/homology/metaPhOrs/ceOrthologs.tsv")

protid.to.gene = ce.metaphors.gene$gene.name
names(protid.to.gene) = ce.metaphors.gene$protid


ce.ortho = read.table(
  "/home/jburdick/data/ftp/phylomedb.org/metaphors/release-201405/orthologs/6239.txt.gz",
  sep="\t", comment="", header=TRUE, as.is=TRUE)
names(ce.ortho)[1] = "taxid1"

# Gets which genes have a Ce ortholog in another species.
# Args:
#   taxid2 - the numerical ID of the other species
# Returns: vector, indexed by Ce gene, of most similar
#   ortholog's score (in [0, 1]; 0 if there is no ortholog)
get.ortho = function(taxid2) {
  s = rep(0, nrow(ce.metaphors.gene))
  names(s) = ce.metaphors.gene[ , "gene.name" ]
  ce.ortho.1 = ce.ortho[ ce.ortho$taxid2 == taxid2 , ]

  ce.ortho.1 = ce.ortho.1[ ce.ortho.1$protid1 %in% names(protid.to.gene) , ]

  ce.ortho.1$gene.name = protid.to.gene[ ce.ortho.1$protid1 ]
  s[ ce.ortho.1$gene.name ] = ce.ortho.1$CS

  s
}


r = get.ortho(7227)







