# Summary of expression for each cell.

source("git/utils.r")
source("git/data/worm/embryodb.r")

movie.list = read.table("git/unmix/image/sort/movieList.tsv",
  sep="\t", header=TRUE, as.is=TRUE)

# XXX this should probably be in git
load("git/unmix/unmix_comp/data/tree_utils.Rdata")

# Sets each cell to the maximum along the lineage path.
lineage.max = function(x) {
  # set each cell to max. of itself and parent
  lt = lin.triples
  x1 = x
  x1[ lt$daughter.1 ] =
    pmax(x1[ lt$parent ], x1[ lt$daughter.1 ])
  x1[ lt$daughter.2 ] =
    pmax(x1[ lt$parent ], x1[ lt$daughter.2 ])

  # if there was no change, we're done
  if (all(x1 == x))
    x1
  else
    # otherwise, make recursive call
    lineage.max(x1)
}

# Summary of per-cell expression.
# Args:
#   series.name - name of the series
# Returns: intensity of that marker for each gene
per.cell.summary = function(series.name) {
  scd = read.embryodb.dat.file(series.name)

  per.cell = c(by(scd$blot, scd$cell,
    function(x) as.numeric(quantile(x, 0.9))))
  per.cell[ "P0" ] = 0
  per.cell = per.cell[ lin.node.names ]

  per.cell
}

# Averages the numbers for each gene.
per.gene.summary = function(gene) {
  x = 0
  series = movie.list[ movie.list$gene == gene, "series" ]
  n = length(series) 
  for(s in series) {
    x = x + per.cell.summary(s)
  }
  x / n
}

r = NULL
for(g in unique(movie.list$gene)) {
  write.status(g)
  r = rbind(r, g = per.gene.summary(g))
}
rownames(r) = unique(movie.list$gene)
rownames(r) = sub("ceh-26", "pros-1", rownames(r))
r = r[ sort(rownames(r)) , ]

write.tsv(r, "git/unmix/image/exprPerCell.tsv")

