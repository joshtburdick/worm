# Adds various sorts of annotation to the TreeView files.

source("git/utils.r")

source("git/data/worm/available_strains.r")

# Get list of genes which have TransgenOme clones
transgenome.genes = {
  tg.table = read.table(
    gzfile("data/genetic/TransgenOmeGenes.tsv.gz"),
    sep="\t", header=TRUE, as.is=TRUE)
  unique(c(tg.table[,5], tg.table[,6]))
}
tg = rep("TG", length(transgenome.genes))
names(tg) = transgenome.genes

# Gets short descriptions of genes.
func.descr = {
  fd0 = read.table(gzfile(
    "data/wormbase/c_elegans.PRJNA13758.WS240.functional_descriptions.txt.gz"),
    sep="\t", skip=4, quote="", fill=TRUE, as.is=TRUE)
  fd0 = fd0[,c(2,7)]
  fd0 = fd0[ !duplicated(fd0[,1]) , ]
  fd1 = fd0[,2]
  names(fd1) = fd0[,1]
  fd1[ fd1 == "not known" ] = NA
  fd1
}

# Anatomy annotation
wb.anatomy = read.table(
  gzfile("data/wormbase/anatomyTerm.ws220.tsv.gz"),
  header=TRUE, sep="\t", quote="", as.is=TRUE)
wb.anatomy = wb.anatomy[ wb.anatomy$Anatomy.Term != "" , ]
anatomy.by.gene = c(by(wb.anatomy$Anatomy.Term,
  wb.anatomy$Gene.Public.Name,
  function(x) paste(sort(unique(as.vector(x))), collapse=",")))
anatomy.by.gene = sapply(anatomy.by.gene,
  function(x) if (nchar(x) >= 140) paste0(substr(x,1,140), "...") else x)

# Adds annotation to a .cdt file.
# Args:
#   input.cdt, output.cdt - .cdt files to read and write
#   annotation - the annotation, as a list of named vectors
# Side effects: writes out output.cdt
add.annotation = function(input.cdt, output.cdt, annotation) {

  r = read.table(input.cdt, sep="\t", as.is=TRUE)
  # get the gene names
  g = r[,2]

  a = matrix(nrow=length(g), ncol=length(annotation))
  rownames(a) = c(1,2,3,g[-c(1:3)])
  colnames(a) = names(annotation)
  # put names of annotation at top
  a[1,] = names(annotation)
  for(n in names(annotation)) {
    g1 = intersect(g, rownames(a))
    a[g1,n] = annotation[[n]][g1]
  }

  r1 = cbind(r[,c(1:3)], a, r[,-c(1:3)])

  write.table(r1, output.cdt,
    sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, na="")
}

# Version of this which works with clustering output,
# as currently defined.
add.annotation.to.clustering = function(cluster.dir) {

  # read in clustering
  cl1 = read.tsv(paste0(cluster.dir, "/clusters.tsv"))
  colnames(cl1) = c("gene", "set")
  cl = cl1$set
  names(cl) = cl1$gene

  annotation = list("Functional description" = func.descr,
    "Anatomy terms" = anatomy.by.gene,
    "cluster" = cl,  # naming this "GROUP" seems to remove dendrogram
    "Strain" = wb.transgene.by.gene,
    "TransGenome" = tg)

  add.annotation(paste0(cluster.dir, "/clusterUnannotated.cdt"),
    paste0(cluster.dir, "/cluster.cdt"),
    annotation)

}

add.annotation.to.clustering("git/cluster/hierarchical/hier.300.clusters/")



