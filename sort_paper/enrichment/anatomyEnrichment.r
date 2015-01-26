# Computes how much anatomy terms (and other known clusterings)
# are enriched in sort fractions and clusters.

source("git/utils.r")
source("git/data/name_convert.r")
source("git/sort_paper/enrichment/hyperg.r")

source("git/sort_paper/FACS/enrichedInFraction.r")

load("git/data/wormbase/anatomy.ontology.group.Rdata")

wb.anatomy = read.table(
  gzfile("data/wormbase/anatomyTerm.ws220.tsv.gz"),
  header=TRUE, sep="\t", quote="", as.is=TRUE)

wb.anatomy.ontology = read.table(
  gzfile("data/wormbase/anatomyOntology.ws220.tsv.gz"),
  header=TRUE, sep="\t", quote="", as.is=TRUE)

anatomy.term.to.name =
  by(wb.anatomy$Anatomy.Term, wb.anatomy$Anatomy.Term.ID,
    function(x) as.vector(x)[1])

if (FALSE) {
# adapt anatomy ontology to be what hyperg.test.groups expects
ao = data.frame(gene=rename.gene.name.vector(wb.anatomy$Gene.Public.Name),
  group=wb.anatomy$Anatomy.Term.ID,
  group.name=wb.anatomy$Anatomy.Term, stringsAsFactors=FALSE)
}

# information about expression clusters from WormBase
load("git/data/wormbase/expr.cluster.Rdata")
wb.cluster = expr.cluster
wb.cluster$gene = rename.gene.name.vector(wb.cluster$gene)
wb.cluster = unique(wb.cluster)
wb.cluster$group.name = expr.cluster.descr[ wb.cluster$group ]

# gene.groups = rbind(ao, wb.cluster)
gene.groups = rbind(ao.group, wb.cluster)

# XXX number of genes affects the significance here
# this number is from "git/sort_paper/plot/numEnrichedInFractions.r"
num.genes = 15683

system("mkdir -p git/sort_paper/enrichment/anatomyEnrichment/")

# compute what's enriched in clusters (old, deprecated version)
if (FALSE) {
for (f in list.files("git/cluster/hierarchical/")) {
  cat(f, "\n")
  cl1 = read.tsv(paste0("git/cluster/hierarchical/", f, "/clusters.tsv"))
  colnames(cl1) = c("gene", "set")
  enrich = hyperg.test.groups.many(gene.groups, cl1, num.genes)
  enrich = enrich[ enrich$p.corr <= 0.05 , ]
  write.tsv(enrich, paste0("git/sort_paper/enrichment/anatomyEnrichment/",
    f, ".tsv"))
}
}

# compute what's enriched in clusters
compute.cluster.enrichment = function() {

  system("mkdir -p git/sort_paper/enrichment/wormbaseCluster/")
  system("mkdir -p git/sort_paper/enrichment/anatomyEnrichment/")

  for (f in list.files("git/cluster/hierarchical/")) {
    cat(f, "\n")
    cl1 = read.tsv(paste0("git/cluster/hierarchical/", f, "/clusters.tsv"))
    colnames(cl1) = c("gene", "set")

    r = hyperg.test.groups.many.faster(unique(wb.cluster), cl1, num.genes)
    save(r, file=paste0("git/sort_paper/enrichment/wormbaseCluster/",
      f, ".Rdata"))

    r = hyperg.test.groups.many.faster(unique(ao.group), cl1, num.genes)
    save(r, file=paste0("git/sort_paper/enrichment/anatomyEnrichment/",
      f, ".Rdata"))
  }
}

# compute what's enriched in sort fractions (old version)
if (FALSE) {

  # compute what's enriched in sort fractions
  source("git/sort_paper/plot/numEnrichedInFractions.r")

  sort.fractions = grep("(\\+|singlets)", colnames(r.sort.only.averaged), invert=TRUE, value=TRUE)
  r = NULL
  for(s in sort.fractions) {

    r1 = as.matrix(r.sort.only.averaged)[ , s ]
    # "background" genes, neither enriched nor depleted
    g.background = names(r1)[ abs(r1) < 2 ]
    g.enriched = names(r1)[ r1 >= 2 ]
    g.depleted = names(r1)[ r1 <= 2 ]

    cl1 = data.frame(gene = c(g.enriched, g.background),
      set = c(rep(paste(s, "enriched"), length(g.enriched)),
        rep("background", length(g.background))))
    r1 = hyperg.test.groups.many(gene.groups, cl1, num.genes)
    r = rbind(r, r1)

    cl1 = data.frame(gene = c(g.depleted, g.background),
      set = c(rep(paste(s, "depleted"), length(g.depleted)),
        rep("background", length(g.background))))
    r1 = hyperg.test.groups.many(gene.groups, cl1, num.genes)
    r = rbind(r, r1)

  }
  r = r[ r$set != "background" , ]
  r$p[ is.nan(r$p) ] = 1        # ??? is this legit?
  r$p.corr = p.adjust(r$p, method="fdr")
  r = r[ r$p.corr <= 0.05 , ]

  write.tsv(r, paste0("git/sort_paper/enrichment/anatomyEnrichment/",
    "sortFraction.tsv"))
}

# enrichment of genes in sort fractions
sort.fraction.anatomy.enrichment = function() {
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

  r = hyperg.test.groups.many.faster(unique(ao.group), cl1, num.genes)
  save(r, file=paste0("git/sort_paper/enrichment/anatomyEnrichment/",
    "facs", ".Rdata"))

  r = hyperg.test.groups.many.faster(unique(wb.cluster), cl1, num.genes)
  save(r, file=paste0("git/sort_paper/enrichment/wormbaseCluster/",
    "facs", ".Rdata"))
}

# enrichment of genes in lineages
if (FALSE) {
  le = read.tsv("git/sort_paper/unmix/lineageEnriched.tsv")
  r = hyperg.test.groups.many(gene.groups, le, num.genes)
  r$p[ is.nan(r$p) ] = 1
  r$p.corr = p.adjust(r$p, method="fdr")
  r = r[ r$p.corr <= 0.05 , ]
  write.tsv(r, "git/sort_paper/enrichment/anatomyEnrichment/lineageEnrichedAnatomyTerms.tsv")
}

# ... and in tissues
if (FALSE) {
  le = read.tsv("git/sort_paper/unmix/tissueEnriched.tsv")
  r = hyperg.test.groups.many(gene.groups, le, num.genes)
  r$p[ is.nan(r$p) ] = 1
  r$p.corr = p.adjust(r$p, method="fdr")
  r = r[ r$p.corr <= 0.05 , ]
  write.tsv(r, "git/sort_paper/enrichment/anatomyEnrichment/tissueEnrichedAnatomyTerms.tsv")
}

# XXX testing with ao, instead of gene groups
# foo = hyperg.test.groups.many.faster(unique(ao), le, num.genes)
# baz = foo$r
# baz["WBbt:0003991","ABalap_enriched",]

# testing with one clustering
if (FALSE) {
cl1 = read.tsv(paste0("git/cluster/hierarchical/hier.300.clusters/clusters.tsv"))
colnames(cl1) = c("gene", "set")
foo = hyperg.test.groups.many.faster(unique(ao.group), cl1, num.genes)
}

sort.fraction.anatomy.enrichment()

if (TRUE) {
  compute.cluster.enrichment()
}


