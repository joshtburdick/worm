# Finds noncoding genes which are in a cluster with
# anatomy annotation (from an enriched set of genes.)

source("git/utils.r")

# statistics about cluster expression and tissue specificity
cluster.ts.stats = read.tsv(
  "git/sort_paper/cluster/annotation/clusterTissueSpecificity.tsv")

# clusters with non-coding RNAs
cl = read.tsv(paste0("git/cluster/hierarchical/hier.300.clusters/clusters.tsv"))

noncoding.cl = cl[ grep("^(anr|linc)-", rownames(cl)) , ]
colnames(noncoding.cl)[1] = "non-coding gene"

# reads per million (for computing max. expression of noncoding genes)
rpm = read.tsv("git/cluster/readsPerMillion.tsv")
rpm = rpm[ , !grepl("^(HS|RNAi)", colnames(rpm)) ]
max.rpm = apply(rpm, 1, max)

noncoding.cl$max.rpm = max.rpm[ noncoding.cl$"non-coding gene" ]

# clusters with some anatomy annotation enriched
anatomy.annotated = read.tsv(
  "git/sort_paper/enrichment/clustersWithAnatomyAnnotation.tsv")

a = merge(noncoding.cl, anatomy.annotated)

a = a[ a$max.rpm >= 1 , ]
a = a[ order(a$max.rpm, decreasing=TRUE) , ]
a = a[ , c("non-coding gene", "max.rpm", "cluster", "group", "p.corr", "group.name") ]

# some stats on this
cat(paste("num expressed anrRNAs with annotation =",
  sum(a$max.rpm >= 1 & grepl("^anr-", a$"non-coding gene")), "\n"))
cat(paste("num expressed lincRNAs with annotation =",
  sum(a$max.rpm >= 1 & grepl("^linc-", a$"non-coding gene")), "\n"))
cat(paste("num clusters containing non-coding genes =",
  length(table(a$cluster)), "\n"))

# noncoding genes in clusters which seem expressed and tissue-specific
noncoding.unannotated =
  noncoding.cl[ ! (noncoding.cl$cluster %in% a$cluster) , ]
noncoding.unannotated =
  noncoding.unannotated[ noncoding.unannotated$max.rpm >= 1,]
noncoding.unannotated =
  noncoding.unannotated[
    cluster.ts.stats[ noncoding.unannotated$cluster , "expr.tissue.spec" ] , ]

cat("num unannotated:",
  length(unique(noncoding.unannotated$"non-coding gene")),
  "ncRNAs in",
  length(unique(noncoding.unannotated$cluster)),
  "clusters\n")


# slight renaming
for(i in 1:nrow(a)) {
  if (grepl("WBPaper00037950", a[i,"group"])) {
    name1 = gsub("WBPaper00037950:", "", a[i,"group"])
    name1 = gsub("-", " ", name1)
    name1 = gsub("_", " ", name1)
    a[i,"group.name"] = name1
  }
}

# complete table
a = merge(noncoding.cl, anatomy.annotated)

a = a[ a$max.rpm >= 1 , ]
a = a[ order(a$max.rpm, decreasing=TRUE) , ]
a = a[ , c("non-coding gene", "max.rpm", "cluster", "group", "p.corr", "group.name") ]

noncoding.anatomy.annotated = a
write.tsv(a,
  "git/sort_paper/noncoding/noncodingInAnnotatedCluster.tsv")

