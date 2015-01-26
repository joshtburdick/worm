# Anatomy enrichment, for a smaller set of broader
# anatomy terms.

source("git/utils.r")
source("git/data/name_convert.r")
source("git/sort_paper/enrichment/hyperg.r")

wb.anatomy = read.table(
  gzfile("data/wormbase/anatomyTerm.ws220.tsv.gz"),
  header=TRUE, sep="\t", quote="", as.is=TRUE)

anatomy.term.to.name =
  by(wb.anatomy$Anatomy.Term, wb.anatomy$Anatomy.Term.ID,
    function(x) as.vector(x)[1])

# adapt anatomy ontology to be what hyperg.test.groups expects
ao = data.frame(gene=rename.gene.name.vector(wb.anatomy$Gene.Public.Name),
  group=wb.anatomy$Anatomy.Term.ID,
  group.name=wb.anatomy$Anatomy.Term, stringsAsFactors=FALSE)

# subset to a smaller set of terms
num.genes.per.anatomy.term = sort(table(ao$group.name), decreasing=TRUE)
anatomy.term.1 = num.genes.per.anatomy.term[1:15]
ao = ao[ao$group.name %in% names(anatomy.term.1) , ] 

# clobber some anatomy terms with less-specific terms
ao[ao$group.name %in% c("head neuron", "nerve ring", "nervous system",
  "neuron", "tail neuron", "ventral cord neuron"), "group.name"] = "neuron"
ao[ao$group.name %in% c("body wall musculature", "vulval muscle"),
  "group.name"] = "muscle"
ao[ao$group.name %in% c("reproductive system", "spermatheca", "vulva"),
  "group.name"] = "reproductive"

ao = unique(ao)

# possibly restrict to genes only annotated with a few terms
specific.gene = {
  gene.count = table(ao$gene)
  names(gene.count)[ gene.count <= 2 ]
}
ao = ao[ ao$gene %in% specific.gene , ]

# information about clusters
wb.cluster = read.table(gzfile("data/wormbase/gene_expr_cluster.tsv.gz"),
  sep="\t", header=TRUE, as.is=TRUE)
wb.cluster = data.frame(gene = rename.gene.name.vector(wb.cluster[,2]),
  group = wb.cluster[,4], group.name = wb.cluster[,5])

gene.groups = rbind(ao, wb.cluster)

# XXX number of genes affects the significance here
# this number is from "git/sort_paper/plot/numEnrichedInFractions.r"
num.genes = 15683

system("mkdir -p git/sort_paper/enrichment/anatomyEnrichment/")

# ... and in tissues
if (TRUE) {
  e = read.tsv("git/sort_paper/unmix/tissueEnriched.tsv")
  e = e[ grep("enriched", e$set), ]
  r = hyperg.test.groups.many(ao, e, num.genes)
  r$p[ is.nan(r$p) ] = 1
  r$p.corr = p.adjust(r$p, method="fdr")
  r = r[ r$p.corr <= 0.05 , ]
  write.tsv(r, "git/sort_paper/enrichment/anatomyEnrichment/tissueEnrichedAnatomyTermsSmall.tsv")
}

# enrichment of things in lineages
if (TRUE) {
  e = read.tsv("git/sort_paper/unmix/lineageEnriched.tsv")
  lin.enriched = e
  e = e[ grep("enriched", e$set), ]
  r = hyperg.test.groups.many(ao, e, num.genes)
  r$p[ is.nan(r$p) ] = 1
  r$p.corr = p.adjust(r$p, method="fdr")
  r = r[ r$p.corr <= 0.05 , ]
  write.tsv(r, "git/sort_paper/enrichment/anatomyEnrichment/lineageEnrichedAnatomyTermsSmall.tsv")
}

