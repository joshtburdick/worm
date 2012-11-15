# Creates a sort matrix from CD files.
# This doesn't include an estimate of cell volume.

source("R/lineage/embryodb.r")

# load("R/lineage/tree_utils.Rdata")
source("R/lineage/tree_utils.r")

# cd = read.embryodb.dat.file("20111011_L1", file.pattern="CD")

# Gets sort matrix entries for a given movie.
# Args:
#   series.name - name of series to use
#   threshold - intensity threshold
# Returns: one row of the sort matrix.
get.sort.matrix.line = function(series.name, threshold) {

  cd = read.embryodb.dat.file(series.name, file.pattern="CD")

  r = c(by(cd$blot >= threshold, cd$cell, mean))

  extend.expr = function(r) {

    # find cases in which leaf isn't present, but parent is
    a = setdiff(names(parent.of), names(r))
    a = a[ !is.na( r[ parent.of[a] ] ) ]

    # extend expression
    r[ a ] = r[ parent.of[a] ]
    r = 1 * (r[ !is.na(r) ] >= 0.1)

    r
  }

  for(i in 1:10) {
#    print(length(r))
    r = extend.expr(r)
  }

  r
}

# Gets the sort matrix for all of the movies
# (possibly including multiple movies per gene.)
# Args: data frame with columns
#   gene - name of the gene
#   series - name of the series
#   threshold - threshold for intensities
# Returns: list with elements:
#   movie.sort.matrix - matrix with one row per series,
#     one column per gene
#   sort.matrix - same, but with all movies for one gene
#     averaged together
get.sort.matrix = function(movie.list) {
  r = NULL

  for(i in 1:nrow(movie.list)) {
    a = movie.list[i,]
    cat(a$series, "")
    m1 = get.sort.matrix.line(a$series, a$threshold)
    m1 = m1[ lin.node.names ]
    r = rbind(r, m1)
  }
  rownames(r) = movie.list[,"series"]
  r[ is.na(r) ] = 0

  list(movie.sort.matrix = r,
    sort.matrix = apply(r, 2,
      function(x) c(by(x, movie.list$gene, mean))))
}

movie.list = read.table("git/unmix/image/sort/movieList.tsv",
  sep="\t", header=TRUE, as.is=TRUE)

r = get.sort.matrix(movie.list)

print(round(apply(r$sort.matrix, 1, sum), 1))

write.table(r$sort.matrix, file="git/unmix/image/sort/sortMatrix.tsv",
  sep="\t", row.names=TRUE, col.names=NA)

