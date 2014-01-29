# Comparing known groups of coregulated genes with clustering.

source("git/utils.r")

# FIXME: construct a matrix with indices of each cluster,
# from all of the clusterings?

# for now, just one of them
cl = {
  cl1 = read.table("git/cluster/hierarchical/hier.200.clusters/clusters.tsv",
    sep="\t", header=TRUE, as.is=TRUE)
  cl = cl1[,3]
  names(cl) = cl1[,2]
  cl
}

# cilia-related genes, from Swoboda 2012
cilia.genes =
  read.table("data/expression/Swoboda2012_RFX_targets.tsv",
    sep="\t", header=TRUE, as.is=TRUE)[,"name"]

# dopaminergic target genes, from Doitsidou et al 2013
dope.target.genes =
  c("dat-1", "cat-4", "cat-2", "bas-1", "asic-1")

# pharyngeal gland cell genes, from Ghai et al 2012
pha.gland.genes = rownames(read.tsv(
  "data/expression/gaudet/Ghai2012_pharyngeal_gland_genes.tsv"))

cat("\nclusters of", length(cilia.genes), "cilia genes:\n")
print(table(cl[cilia.genes]))

cat("\nclusters of", length(dope.target.genes), "dopaminergic target genes:\n")
print(table(cl[dope.target.genes]))

cat("\nclusters of", length(pha.gland.genes), "pharyngeal gland genes:\n")
print(table(cl[pha.gland.genes]))

