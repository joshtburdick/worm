# Utilities for computing hypergeometric enrichment.

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

# was:
#  p = phyper(length(g) - num.in.group,
#    num.genes - length(group), length(group), length(g),
#    lower.tail=TRUE)
# new version, based on
# http://stackoverflow.com/questions/8382806/r-hypergeometric-test-phyper
  p = 1-phyper(num.in.group-1, length(group),
    num.genes-length(group), length(g))

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
      num.in.group = p$num.in.group, p.corr=NA, genes = p$genes))
  }

  r$p.corr = p.adjust(r$p, method="BH")
  r
}

