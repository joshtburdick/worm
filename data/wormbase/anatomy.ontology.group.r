# Groups together all descendents of particular anatomy
# ontology terms.

wb.anatomy = read.table(
  gzfile("data/wormbase/anatomyTerm.ws220.tsv.gz"),
  header=TRUE, sep="\t", quote="", as.is=TRUE)

wb.anatomy.ontology = read.table(
  gzfile("data/wormbase/anatomyOntology.ws220.tsv.gz"),
  header=TRUE, sep="\t", quote="", as.is=TRUE)



anatomy.term.parent = strsplit(wb.anatomy.ontology[,3], " \\| ")
names(anatomy.term.parent) = wb.anatomy.ontology[,1]
anatomy.term.parent =
  anatomy.term.parent[ sapply(anatomy.term.parent, length) > 0 ]

if (FALSE) {
ao.table = NULL
for(a in names(anatomy.term.parent))
  ao.table = rbind(ao.table,
    data.frame(a, anatomy.term.parent[[a]], stringsAsFactors=FALSE))
}

# anatomy.term.child =
#   by(ao.table[,1], ao.table[,2], function(x) unique(as.character(x)))

# a = relation(graph=ao.table)


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

# Computes a transitive closure. This is essentially a sparse
# version of Floyd-Warshall. (Deprecated: using "relation" instead.)
# (No, maybe not.)
transitive.closure.list = function(r) {

  print(object.size(r))

  # function which adds things onto a list
  closure.step = function(x) {
    a = unique(c(x, r[[x]], sapply(r[[x]], function(y) r[[y]])))
    a = a[ !is.null(a) ]
    a
  }

  r1 = sapply(names(r), closure.step)

  # are we done?
  if (length(c(r1, recursive=TRUE)) == length(c(r, recursive=TRUE)))
    r
  else
    # no; make recursive call
    transitive.closure.list(r1)
}

# foo = by(c(1:10), c(1,1,1,2,2,2,2,3,3,3), c, simplify=FALSE)

# foo = transitive.closure.list(anatomy.term.child)

a = relation.to.matrix(anatomy.term.parent)
a1 = refl.trans.closure.matrix(a)

