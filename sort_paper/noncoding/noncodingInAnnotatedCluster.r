# Finds noncoding genes which are in a cluster with
# anatomy annotation (from an enriched set of genes.)

source("git/utils.r")

# clusters with non-coding RNAs
cl = read.tsv(paste0("git/cluster/hierarchical/hier.300.clusters/clusters.tsv"))

noncoding.cl = cl[ grep("^(anr|linc)-", rownames(cl)) , ]
colnames(noncoding.cl)[1] = "non-coding gene"

# clusters with some anatomy annotation enriched
anatomy.annotated = read.tsv(
  "git/sort_paper/enrichment/clustersWithAnatomyAnnotation.tsv")

noncoding.anatomy.annotated = merge(noncoding.cl, anatomy.annotated)

# reads per million (for computing max. expression)
rpm = read.tsv("git/cluster/readsPerMillion.tsv")
rpm = rpm[ , !grepl("^(HS|RNAi)", colnames(rpm)) ]
max.rpm = apply(rpm, 1, max)

noncoding.anatomy.annotated$max.rpm =
  max.rpm[ noncoding.anatomy.annotated$"non-coding gene" ]
noncoding.anatomy.annotated = noncoding.anatomy.annotated[
  order(noncoding.anatomy.annotated$max.rpm, decreasing=TRUE) , ]
noncoding.anatomy.annotated = noncoding.anatomy.annotated[ ,
  c("non-coding gene", "max.rpm", "cluster", "group", "p.corr", "group.name") ]

write.tsv(noncoding.anatomy.annotated,
  "git/sort_paper/noncoding/noncodingInAnnotatedCluster.tsv")

