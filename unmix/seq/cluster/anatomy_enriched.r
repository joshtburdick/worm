# Looks for anatomy terms which are overrepresented in one fraction or another.

wb.anatomy = read.table(
  gzfile("data/wormbase/anatomyTerm.ws220.tsv.gz"),
  header=TRUE, sep="\t", quote="", as.is=TRUE)

anatomy.term.to.name = by(wb.anatomy$Anatomy.Term, wb.anatomy$Anatomy.Term.ID,
  function(x) as.vector(x)[1])

# r = as.matrix(read.table("git/unmix/seq/quant/readsPerMillion.tsv",
#   header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))
r = as.matrix(read.table("git/unmix/seq/cluster/readsNormalized.tsv",
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

# lr = log2(1 + r)
# control = (lr[,"cnd-1_ungated"] + lr[,"pha-4_ungated"]) / 2
# control = (lr[,"cnd-1_singlets"] + lr[,"pha-4_singlets"]) / 2
# le = lr - control
le = r[,1:25]

# Hypergeometric test for enrichment.
# Args:
#   g - set of genes
#   group - a group of genes
#   num.genes - the total number of genes
# Returns: p-value for enrichment of genes in group
#   in the set g.
hyperg.test.enrichment = function(g, group, num.genes) {
  # make sure groups of genes are sets
  g = unique(g)
  group = unique(group)

  # count genes in group
  num.in.group = length(intersect(g, group))

  p = phyper(length(g) - num.in.group,
    num.genes - length(group), length(group), length(g),
    lower.tail=TRUE)
  list(p = p, num.intersect = length(intersect(g, group)),
    num.in.group = length(group),
    genes = if (length(intersect(g, group)) > 0) paste(intersect(g,group), collapse=" ") else "")
}

# Tests for enrichment of a list of genes.
# Args:
#   groups - a data frame with columns
#     "gene", "group", and "group.name"
#   gene.list - a list of genes
#   num.genes - the assumed total number of genes
# Returns: a data frame of results
hyperg.test.groups = function(groups, gene.list, num.genes) {
  r = NULL

  for(gr in unique(groups$group)) {
# cat(gr, "")
    g1 = groups[ groups$group==gr, ]
    p = hyperg.test.enrichment(gene.list, g1$gene, num.genes)
    r = rbind(r, data.frame(group=gr, group.name=g1$group.name[1],
      p = p$p, num.intersect = p$num.intersect,
      num.in.group = p$num.in.group, genes = p$genes))
  }

  r$p.corr = p.adjust(r$p, method="BH")
  r
}

# computes anatomy terms enriched in each fraction
compute.anatomy.enriched = function() {
  system("mkdir -p git/unmix/seq/cluster/anatomy_enriched/")
  for(g in colnames(lr)) {
    cat(g, "")
    r = enriched.t.test(lr[,g], control)
    r = r[ r$p.corr <= 0.05 & r$statistic > 0, ]
    r = r[ order(r$p.corr), ]
    write.table(r, sep="\t", col.names=NA,
      file=paste("git/unmix/seq/cluster/anatomy_enriched/", g, ".tsv", sep=""))
  }
}

# similar comparison for, e.g., negatives
# Args:
#   a, b - names of two samples
# Side effects: writes comparisons
compute.anatomy.enriched.1 = function(a, b) {
  system("mkdir -p git/unmix/seq/cluster/anatomy_enriched/")
  cat(a, b, "\n")

  r = enriched.t.test(lr[,a], lr[,b])
  r = r[ r$p.corr <= 0.05 & r$statistic > 0, ]
  r = r[ order(r$p.corr), ]
  name = paste(a, "rel_to", b, sep="_")
  write.table(r, sep="\t", col.names=NA,
    file=paste("git/unmix/seq/cluster/anatomy_enriched/", name, ".tsv", sep=""))
}

# Computes anatomy enrichment for the negatives.
compute.anatomy.enriched.negatives = function() {
  compute.anatomy.enriched.1("ceh-36p", "ceh-36m")

  compute.anatomy.enriched.1("cnd-1p8_19", "cnd-1m")
  compute.anatomy.enriched.1("cnd-1_12_14", "cnd-1m")
  compute.anatomy.enriched.1("cnd-1p1_4", "cnd-1m")

  compute.anatomy.enriched.1("pha-4p9_1", "pha-4m")
  compute.anatomy.enriched.1("pha-4p12_9", "pha-4m")

  compute.anatomy.enriched.1("cnd-1_ungated", "pha-4_ungated")
  compute.anatomy.enriched.1("cnd-1_singlets", "pha-4_singlets")
}

# adapt anatomy ontology to be what hyperg.test.groups expects
ao = data.frame(gene=wb.anatomy$Gene.Public.Name,
  group=wb.anatomy$Anatomy.Term.ID,
  group.name=wb.anatomy$Anatomy.Term, stringsAsFactors=FALSE)
num.genes = dim(r)[1]

compute.anatomy.enriched.hypergeometric = function() {
  system("mkdir -p git/unmix/seq/cluster/anatomy_enriched_hyperg/")
#  for(fr in c("ceh-26", "ceh-27", "ceh-36", "ceh-6", "cnd-1",
#    "F21D5.9", "hlh-16", "irx-1", "mir-57", "mls-2", "pal-1",
#    "pha-4", "ttx-3", "unc-130")) {
  files = list.files("git/unmix/seq/cluster/enriched_genes/", pattern=".tsv")
  samples = sub(".tsv", "", files)

  for(fr in samples) {
cat(fr, "")
    try({
    g = read.table(
      paste("git/unmix/seq/cluster/enriched_genes/", fr, ".tsv", sep=""),
      as.is=TRUE, header=FALSE)[,1]
    a = hyperg.test.groups(ao, g, num.genes)
    a = a[a$p < 0.05,]
    a = a[order(a$p),]
    write.table(a, sep="\t", col.names=NA,
      file=paste("git/unmix/seq/cluster/anatomy_enriched_hyperg/", fr, ".tsv", sep=""))
    })
  }
}


# compute.anatomy.enriched()
# compute.anatomy.enriched.negatives()

# compute.anatomy.enriched.hypergeometric()

