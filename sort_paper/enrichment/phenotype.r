# Computes how many phenotypes are enriched in
# sort fractions and clusters.

source("git/utils.r")
source("git/data/name_convert.r")
source("git/sort_paper/enrichment/hyperg.r")
source("git/data/wormbase/refl.trans.closure.r")

source("git/sort_paper/FACS/enrichedInFraction.r")

output.dir = "git/sort_paper/enrichment/phenotype/allResults/"
system(paste0("mkdir -p ", output.dir))

# XXX number of genes affects the significance here
# this number is from "git/sort_paper/plot/numEnrichedInFractions.r"
num.genes = 15683

# group phenotypes according to ontology
# this already includes the reflexive, transitive closure
phenotype.descendent = read.table(
  gzfile("data/wormbase/phenotype_ontology_WS220.tsv.gz"),
  sep = "\t", header=TRUE, as.is=TRUE)
colnames(phenotype.descendent) =
  c("Parent.Phenotype.ID", "Phenotype.ID")
# pheno.m = edges.to.matrix(
#   phenotype.descendent[,2], phenotype.descendent[,1])
# pheno.ontology = refl.trans.closure.matrix(pheno.m)

# information about phenotypes
phenotype = read.tsv(gzfile("data/wormbase/phenotype_WS220.tsv.gz"))

# which gene has which phenotype
phenotype.gene = read.table(
  "data/wormmine/phenotype_gene.tsv.gz",
  sep="\t", header=TRUE, as.is=TRUE)
# use sequence name, if gene name isn't defined
i = phenotype.gene[,2] == ""
phenotype.gene[i,2] = phenotype.gene[i,1]
phenotype.gene = phenotype.gene[,c(3,2)]
colnames(phenotype.gene) = c("Phenotype.ID", "gene")

phenotype.gene = merge(phenotype.descendent, phenotype.gene)
phenotype.gene$phenotype.name =
  phenotype[phenotype.gene$Parent.Phenotype.ID, "Phenotype Name"]

# RNAi phenotypes
rnai = read.tsv("git/data/wormbase/rnai.phenotype.tsv.gz")
colnames(rnai) = c("gene", "Phenotype.ID")
rnai = merge(phenotype.descendent, rnai)
rnai$phenotype.name =
  phenotype[rnai$Parent.Phenotype.ID, "Phenotype Name"]

# combine these
pheno = rbind(
  data.frame(
    gene = phenotype.gene$gene,
    group = phenotype.gene$Parent.Phenotype.ID,
    group.name = phenotype.gene$phenotype.name, stringsAsFactors=FALSE),
  data.frame(gene = rnai$gene,
    group = paste0("RNAi_", rnai$Parent.Phenotype.ID),
    group.name = paste(rnai$phenotype.name, "(RNAi)"), stringsAsFactors=FALSE))

pheno$gene = rename.gene.name.vector(pheno$gene)
pheno = pheno[ pheno$gene != "" , ]
pheno = pheno[ pheno$group != "" , ]
pheno = unique(pheno)


# compute what's enriched in clusters
compute.cluster.enrichment = function() {
  for (f in c("hier.300.clusters", list.files("git/cluster/hierarchical/"))) {
    cat(f, "\n")
    cl1 = read.tsv(paste0("git/cluster/hierarchical/", f, "/clusters.tsv"))
    colnames(cl1) = c("gene", "set")

    r = hyperg.test.groups.many.faster(unique(pheno), cl1, num.genes)
    r[,,"p.corr"][ is.na(r[,,"p.corr"]) ] = 1
    save(r, file=paste0(output.dir, f, ".Rdata"))
  }
}

# enrichment of genes in sort fractions
sort.fraction.enrichment = function() {
  cl1 = NULL

  for(g in names(facs.enriched.depleted)) {
    a = facs.enriched.depleted[[g]]
    if (sum(a) > 0)
      cl1 = rbind(cl1, data.frame(gene = names(a)[a], set = g,
        stringsAsFactors=FALSE))
  }
  cl1$gene = as.character(cl1$gene)
  cl1$set = as.character(cl1$set)

  colnames(cl1) = c("gene", "set")

  r = hyperg.test.groups.many.faster(pheno, cl1, num.genes)
  r[,,"p.corr"][ is.na(r[,,"p.corr"]) ] = 1
  save(r, file=paste0(output.dir, "facs", ".Rdata"))
}

sort.fraction.enrichment()
compute.cluster.enrichment()


