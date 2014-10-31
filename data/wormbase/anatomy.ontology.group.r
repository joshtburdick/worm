# Groups together all descendents of particular anatomy
# ontology terms.

library(Matrix)

source("git/utils.r")
source("git/data/name_convert.r")

wb.anatomy = read.table(
  gzfile("data/wormbase/anatomyTerm.ws220.tsv.gz"),
  header=TRUE, sep="\t", quote="", as.is=TRUE)

wb.anatomy.ontology = read.table(
  gzfile("data/wormbase/anatomyOntology.ws220.tsv.gz"),
  header=TRUE, sep="\t", quote="", as.is=TRUE)

# anatomy.term.to.name =
#   by(wb.anatomy$Anatomy.Term, wb.anatomy$Anatomy.Term.ID,
#     function(x) as.vector(x)[1])
anatomy.term.to.name = wb.anatomy.ontology$Anatomy.Term
names(anatomy.term.to.name) = wb.anatomy.ontology$Anatomy.Term.ID
anatomy.term.to.name[ anatomy.term.to.name == "" ] =
  wb.anatomy.ontology[ anatomy.term.to.name == "", "Anatomy.Term.ID" ]

anatomy.term.parent = strsplit(wb.anatomy.ontology[,3], " \\| ")
names(anatomy.term.parent) = wb.anatomy.ontology[,1]
anatomy.term.parent =
  anatomy.term.parent[ sapply(anatomy.term.parent, length) > 0 ]

# adapt anatomy ontology to be what hyperg.test.groups expects
ao = data.frame(gene=rename.gene.name.vector(wb.anatomy$Gene.Public.Name),
  group=wb.anatomy$Anatomy.Term.ID,
  group.name=wb.anatomy$Anatomy.Term, stringsAsFactors=FALSE)

# Converts a relation (in adjacency-list form) to a logical matrix.
relation.to.matrix = function(r) {
  a = unique(c(names(r), r, recursive=TRUE))
  m = Matrix(FALSE, nrow=length(a), ncol=length(a),
    dimnames = list(a, a))

  for(i in names(r))
    m[ i, r[[i]] ] = TRUE

  m
}

# Computes the reflexive, transitive closure of a relation.
# (Basically a dense version of Floyd-Warshall.)
refl.trans.closure.matrix = function(r) {
  diag(r) = TRUE

  r1 = (r + r %*% r) > 0
  while (any(r1 != r)) {
cat(sum(r), " ")
    r = r1
    r1 = (r + r %*% r) > 0
  }

  r
}

a = relation.to.matrix(anatomy.term.parent)
a1 = refl.trans.closure.matrix(a)

# compute list of genes in each anatomy ontology group
ao.group.l = vector("list", ncol(a1))
for(j in 1:ncol(a1)) {
  g = colnames(a1)[j]

  # human-readable name for this term (if available)
  group.name = anatomy.term.to.name[ g ]
  if (is.na(group.name))
    group.name = g

  write.status(paste(j, g, group.name))
  g1 = rownames(a1)[ a1[,j] ]
  genes = unique( ao[ ao$group %in% g1, "gene" ] )
  if (length(genes) > 0) {
    ag = data.frame(gene = genes, group = g, group.name = group.name,
      row.names = genes, stringsAsFactors = FALSE)
    ao.group.l[[j]] = ag     # ao.group = rbind(ao.group, ag)
  }
}
ao.group = do.call("rbind", ao.group.l)
# ao.group$group.name = anatomy.term.to.name[ ao.group$group ]

ao.term.to.name = anatomy.term.to.name
ao.subterm = a1
save(ao, ao.group, anatomy.term.to.name, ao.subterm,
  file="git/data/wormbase/anatomy.ontology.group.Rdata")

