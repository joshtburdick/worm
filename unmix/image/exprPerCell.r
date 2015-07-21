# Summary of expression for each cell.

source("git/utils.r")
source("git/data/worm/embryodb.r")

movie.list = read.table("git/unmix/image/sort/movieList.tsv",
  sep="\t", header=TRUE, as.is=TRUE)
# arbitrarily picking the first movie per gene
movie.list = movie.list[ ! duplicated(movie.list$gene) , ]

# XXX this should probably be in git
load("git/unmix/unmix_comp/data/tree_utils.Rdata")

# anatomy annotation of these
tissues.per.cell = read.table("data/worm/TissuesPerCell.tsv",
  sep="\t", quote="", header=TRUE, row.names=1, as.is=TRUE)

# Summary of per-cell expression.
# Args:
#   series.name - name of the series
# Returns: intensity of that marker for each cell
series.summary = function(series.name) {
  sca = read.embryodb.dat.file(series.name, file.pattern="SCA")
  sca = sca[sca$cell %in% lin.node.names , ]
  rownames(sca) = sca$cell
  expr = sca[ lin.node.names, "blot" ]
  names(expr) = lin.node.names
  expr[ "P0" ] = 0
  expr
}

r = NULL
for(i in 1:nrow(movie.list)) {
  r = rbind(r, g = series.summary(movie.list[i,"series"]))
}
rownames(r) = movie.list$gene
rownames(r) = sub("ceh-26", "pros-1", rownames(r))
r = r[ sort(rownames(r)) , ]

tissues.per.cell = tissues.per.cell[ colnames(r) , ]

r = data.frame(cell = tissues.per.cell$Cell,
  tissue = tissues.per.cell$Tissue,
  t(r), check.names = FALSE, stringsAsFactors = FALSE)

r[ is.na(r$cell), "cell" ] = ""
r[ is.na(r$tissue), "tissue" ] = ""

write.tsv(r, "git/unmix/image/exprPerCell.tsv")

