# This is the sort matrix, in a table.

source("git/utils.r")

# the sort matrix (unnormalized)
source("git/sort_paper/unmix/sortMatrix.r")
m1 = 1 * (m.unnormalized > 0)

m1 = m1[ rownames(m1) != "ceh-6 (-) hlh-16 (-)" , ]

# anatomy annotation of these
tissues.per.cell = read.table("data/worm/TissuesPerCell.tsv",
  sep="\t", quote="", header=TRUE, row.names=1, as.is=TRUE)
tissues.per.cell = tissues.per.cell[ colnames(m1) , ]

m1 = t(rbind(cell = tissues.per.cell$Cell,
  tissue = tissues.per.cell$Tissue,
  m1))

write.tsv(m1, file="git/sort_paper/unmix/sortMatrix/table.tsv")

