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

a = merge(noncoding.cl, anatomy.annotated)

# reads per million (for computing max. expression)
rpm = read.tsv("git/cluster/readsPerMillion.tsv")
rpm = rpm[ , !grepl("^(HS|RNAi)", colnames(rpm)) ]
max.rpm = apply(rpm, 1, max)

a$max.rpm = max.rpm[ a$"non-coding gene" ]
a = a[ order(a$max.rpm, decreasing=TRUE) , ]
a = a[ , c("non-coding gene", "max.rpm", "cluster", "group", "p.corr", "group.name") ]

# some stats on this
cat(paste("num expressed anrRNAs with annotation =",
  sum(a$max.rpm >= 1 & grepl("^anr-", a$"non-coding gene")), "\n"))
cat(paste("num expressed lincRNAs with annotation =",
  sum(a$max.rpm >= 1 & grepl("^linc-", a$"non-coding gene")), "\n"))

write.tsv(noncoding.anatomy.annotated,
  "git/sort_paper/noncoding/noncodingInAnnotatedCluster.tsv")


