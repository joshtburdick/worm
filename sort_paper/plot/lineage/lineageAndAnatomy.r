# Plots some lineages, together with the anatomy
# terms for those cells.

source("git/utils.r")
source("git/plot/utils.r")
source("git/plot/plot_expr.r")

# the sort matrix (unnormalized)
source("git/sort_paper/unmix/sortMatrix.r")

# annotation of which tissues are in which cell, and coloring
# of them
tissues.per.cell = read.table("data/worm/TissuesPerCell.tsv",
  sep="\t", quote="", header=TRUE, row.names=1, as.is=TRUE)
tissues.per.cell.1 = sapply(lin.node.names,
  function(cell) {
    leaves = names(which(leaf.lineage.matrix[ cell , ] == 1))
    tissues = tissues.per.cell[ leaves, "Tissue" ]
    if (mean(tissues == tissues[1]) == 1)
      tissues[1]
    else
      ""
  })
# various tweaks and simplifications
tissues.per.cell.1[ tissues.per.cell.1 == "Arcade" ] = "Pharynx"
tissues.per.cell.1[ tissues.per.cell.1 == "Nervous" ] = "Neuron"

tissue.to.color = NULL
tissue.to.color["Death"] = "pink"
tissue.to.color["Epidermis"] = "lightblue"
tissue.to.color["Glia"] = "orange"
tissue.to.color["Intestine"] = "darkgreen"
tissue.to.color["Muscle"] = "magenta"
tissue.to.color["Neuron"] = "blue"
tissue.to.color["Pharynx"] = "green"


# unc-130
a = rgb(m.unnormalized["unc-130",], 0, 0)
names(a) = colnames(m.unnormalized)
a["ABplapapppa"] = rgb(1,0,0)    # XXX hack

lineage.to.plot = "ABpr"

cells = get.node.names(get.lineage.by.name(lin, lineage.to.plot), 100)
cells = colnames(m.unnormalized)[ colnames(m.unnormalized) %in% cells ]
cells = cells[ nchar(cells) == 11 ]       # only include terminal cells
cell.col = number.lineage.columns(get.lineage.by.name(lin, lineage.to.plot))

# XXX in order to put later tissue fates higher
cell.col.1 = cell.col[order(nchar(cell.col))]

pdf("git/sort_paper/plot/lineage/lineageAndAnatomy_unc130.pdf",
  width=8, height=4)

plot(1,1, xlim=range(cell.col), ylim=c(540,0),
  type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")

mtext("Tissue fate", side=1, cex=1.5)
mtext(expression(italic("unc-130") * " expression"),
  side=3, cex=1.5, col="red")

par(new=FALSE)
plot.segments.per.cell(a, main, root=lineage.to.plot, times=c(50,440),
      add=TRUE, lwd=3, yaxt="n")

rect(cell.col.1 - 0.5, 460, cell.col.1 + 0.5, 500,
  border=NA,
  col=tissue.to.color[ tissues.per.cell.1[ names(cell.col.1) ] ])
dev.off()






