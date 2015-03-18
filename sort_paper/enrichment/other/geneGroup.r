# Looking to see if various categories of genes are enriched in
# particular clusters.

source("git/utils.r")

source("git/sort_paper/enrichment/hyperg.r")

# using orthology annotation from MetaPhOrs
ortho = as.matrix(read.tsv(gzfile(
  "git/data/homology/metaPhOrs/has.ortholog.tsv.gz")))

# construct table of what's orthologous in each cluster
has.ortho = NULL
for(s in colnames(ortho)) {
  g = rownames(ortho)[ ortho[ , s ] > 0 ]
  a1 = data.frame(gene = g, group = s, group.name = s,
    stringsAsFactors=FALSE)
  has.ortho = rbind(has.ortho, a1)
}

# add on Pristinchous-specific
o = ortho > 0
g1 = o[,"Pristionchus"] & (! (o[,"fly"] | o[,"human"] | o[,"mouse"]))
g1 = names(g1)[g1]
a1 = data.frame(gene = g1, group = "Pp only", group.name = "Pp only",
  stringsAsFactors=FALSE)
has.ortho = rbind(has.ortho, a1)

has.ortho = unique(has.ortho)

# a clustering
cl1 = read.tsv(paste0("git/cluster/hierarchical/hier.300.clusters/clusters.tsv"))
colnames(cl1) = c("gene", "set")

# XXX number of genes affects the significance here
# this number is from "git/sort_paper/plot/numEnrichedInFractions.r"
num.genes = 15683
r = hyperg.test.groups.many.faster(has.ortho, cl1, num.genes)






