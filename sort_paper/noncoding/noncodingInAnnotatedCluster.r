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

# require at least this much expression of the actual non-coding gene
noncoding.cl$max.rpm = max.rpm[ noncoding.cl$"non-coding gene" ]
noncoding.cl = noncoding.cl[ noncoding.cl$max.rpm >= 1 , ]

# clusters with some anatomy annotation enriched
anatomy.annotated = read.tsv(
  "git/sort_paper/enrichment/clustersWithAnatomyAnnotation.tsv")

# slight renaming
for(i in 1:nrow(anatomy.annotated)) {
  if (grepl("WBPaper00037950", anatomy.annotated[i,"group"])) {
    name1 = gsub("WBPaper00037950:", "", anatomy.annotated[i,"group"])
    name1 = gsub("-", " ", name1)
    name1 = gsub("_", " ", name1)
    name1 = sub(" (enriched|CoreEnriched|SelectivelyEnriched)", "", name1)
    anatomy.annotated[i,"group.name"] = name1
  }
}
anatomy.annotated = anatomy.annotated[ , c("cluster", "group", "group.name", "p.corr") ]

# do an "outer join"
a = merge(noncoding.cl, anatomy.annotated, all.x=TRUE)

# for genes without an annotated cluster, require a minimum cluster
# expression and tissue-specificity
# a = a[ (!is.na(a$group)) || cluster.ts.stats[a$cluster, "expr.tissue.spec"] , ]

# sort the table
a = a[ order(is.na(a$group), a$max.rpm, decreasing=TRUE) , ]
a = a[ , c("non-coding gene", "max.rpm", "cluster", "group", "p.corr", "group.name") ]

# some stats on this
cat("num expressed anrRNAs with annotation =",
  sum((!is.na(a$group)) & grepl("^anr-", a$"non-coding gene")), "\n")
cat("num expressed lincRNAs with annotation =",
  sum((!is.na(a$group)) & grepl("^linc-", a$"non-coding gene")), "\n")
annotated.group = a$group
annotated.group = annotated.group[ !is.na(annotated.group) ]
cat("num annotated clusters containing non-coding genes =",
  length(unique(annotated.group)), "\n")

# noncoding genes in clusters which seem expressed and tissue-specific
noncoding.unannotated = a[ is.na(a$group) , ]
cat("num unannotated:",
  length(unique(noncoding.unannotated$"non-coding gene")),
  "ncRNAs in",
  length(unique(noncoding.unannotated$cluster)),
  "clusters\n")

# write out table
a = a[ order(!is.na(a$group), a$max.rpm, decreasing=TRUE) , ]
a = a[ , c("non-coding gene", "max.rpm", "cluster", "group", "p.corr", "group.name") ]
colnames(a) = c("Non-coding gene", "Max. RPM", "Cluster", "Anatomy ID", "Anatomy p", "Anatomy group")
write.tsv(a,
  "git/sort_paper/noncoding/noncodingInAnnotatedCluster.tsv")

