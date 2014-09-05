# Utilities for computing hypergeometric enrichment.

library(Matrix)

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

  # new version, based on
  # http://stackoverflow.com/questions/8382806/r-hypergeometric-test-phyper
  # p = 1-phyper(num.in.group-1, length(group),
  #    num.genes-length(group), length(g))

  # ...except, for reasons of numerical stability (not
  # computing "1 - 0.999999999...whatever"), using this instead
  p = phyper(num.in.group-1, length(group),
    num.genes-length(group), length(g), lower.tail=FALSE)

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
    g1 = groups[ groups$group==gr, ]
    p = hyperg.test.enrichment(gene.list, g1$gene, num.genes)
    r = rbind(r, data.frame(group=gr, group.name=g1$group.name[1],
      p = p$p, num.intersect = p$num.intersect,
      num.in.group = p$num.in.group, p.corr=NA, genes = p$genes))
  }

  r
}

# Tests for enrichment of any of several lists of genes,
# in a large number of groups of genes.
# Args:
#   groups - a data frame of, e.g., anatomy terms, with
#     columns "gene", "group", and "group.name"
#   gene.set - a data.frame with columns "gene" and "set"
#     (e.g. cluster)
#   num.genes - the assumed total number of genes
# Returns: a data frame of results
# XXX somewhat slow and deprecated
hyperg.test.groups.many = function(groups, gene.set, num.genes) {
  gene.set = gene.set[ !is.na(gene.set$set) , ]

  r = NULL

  for(s in sort(unique(gene.set$set))) {
    write.status(s)
    g = gene.set[ gene.set$set==s, "gene" ]
    r1 = hyperg.test.groups(groups, g, num.genes)
    r = rbind(r, data.frame(set = s, r1))
  }

  r$p.corr = p.adjust(r$p, method="fdr")
  r = r[ order(r$p.corr, r$p) , ]
  r
}

# Converts a table to a logical matrix.
table.to.logical.matrix = function(i, j) {
  R = 1:length(unique(i))
  C = 1:length(unique(j))
  names(R) = sort(unique(i))
  names(C) = sort(unique(j))

  a = Matrix(FALSE, nrow=length(R), ncol=length(C))
  rownames(a) = names(R)
  colnames(a) = names(C)

  a[ cbind(R[i],C[j]) ] = TRUE

  a
}

# Faster version of this, using matrices.
# Args:
#   groups - a data frame of, e.g., anatomy terms, with
#     columns "gene", "group", and "group.name"
#   gene.set - a data.frame with columns "gene" and "set"
#     (e.g. cluster)
#   num.genes - the assumed total number of genes
# Returns: an array of results
hyperg.test.groups.many.faster = function(groups, gene.set, num.genes) {

  # XXX avoid problems with cluster names which are numbers
  gene.set$set = as.character(gene.set$set)

  # mapping of groups to human-readable name (assumes that a given
  # anatomy term always has the same human-readable name)
  # XXX this is "non-normalized" ??? make this an arg.?
  group.to.name =
    by(groups$group.name, groups$group,
    function(x) as.vector(x)[1])

  a = table.to.logical.matrix(groups$group, groups$gene)
  b = table.to.logical.matrix(gene.set$gene, gene.set$set)
  group.sizes = apply(a, 1, sum)
  cluster.sizes = apply(b, 2, sum)

  genes = intersect(colnames(a), rownames(b))
  counts = a[,genes] %*% b[genes,]

  r = array(NA, dim=c(nrow(counts), ncol(counts), 5),
    dimnames = list(group=rownames(counts), cluster=colnames(counts),
    stat=c("num.in.group", "num.in.cluster", "num.intersect", "p", "p.corr")))

  for(group in rownames(counts))
    for(cl in colnames(counts)) {
write.status(paste(group, cl))
      num.in.group = group.sizes[group]
      num.in.cluster = cluster.sizes[cl]
      p = phyper(counts[group,cl]-1, num.in.group,
        num.genes-num.in.group, num.in.cluster, lower.tail=FALSE)
      r[ as.character(group) , as.character(cl) , ] =
        c(num.in.group, num.in.cluster, counts[group,cl], p, NA)
  }

  r[,,"p.corr"] = p.adjust(r[,,"p"], method="fdr")

#  list(a = a, b = b, counts = counts, r = r)
  r
}


