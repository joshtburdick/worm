# Plots a summary of which genes were enriched in which
# lineage.

# library(ggplot2)
# library(ggdendro)

library(plotrix)

source("git/utils.r")
source("git/plot/plot_expr.r")
source("git/sort_paper/enrichment/anatomyEnrichmentSmall.r")

tissues.per.cell = read.table("data/worm/TissuesPerCell.tsv",
  sep="\t", quote="", header=TRUE, row.names=1, as.is=TRUE)

# the cells which are included as leaf nodes
cells.23 = c("ABalaa", "ABalap", "ABalpa", "ABalpp", 
  "ABaraa", "ABarap", "ABarpa", "ABarpp",
  "ABplaa", "ABplap", "ABplpa", "ABplpp",
  "ABpraa", "ABprap", "ABprpa", "ABprpp",
  "MSaa", "MSap", "MSpa", "MSpp", "E", "C", "P3")

# Gets counts of known tissue annotations for a lineage.
# XXX requires a certain amount of munging
known.tissue.count = function(lineage) {
  cells = names(which(cell.lineage.matrix[lineage,]>0))
  tissues = tissues.per.cell[ cells , "Tissue" ]

  tissues = tissues[ ! tissues %in% c("", "Death", "Other") ]
  c(hypdermis = sum(tissues %in% c("Epidermis")),
    intestine = sum(tissues %in% c("Intestine")),
    muscle = sum(tissues %in% c("Muscle")),
    neuron = sum(tissues %in% c("Glia", "Nervous")),
    pharynx = sum(tissues %in% c("Pharynx")),
    reproductive = sum(tissues %in% c("Gonad")),
    "seam cell" = 0)
}

# Plots using ggplot (deprecated).
plot.dendro.gg = function() {
  ddata = dendro_data(lin)

  ### Set up a blank theme
  theme_none <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(colour=NA),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank()
    #axis.ticks.length = element_blank()
  )

  # p.lin <- ggplot(segment(ddata)) + 
  #   geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + 
  #   theme_none + theme(axis.title.x=element_blank())
  p.lin = ggdendrogram(lin) + theme_bw()

  ggsave(p.lin, file="git/sort_paper/plot/lineage/lineageEnriched.png",
    width=9, height=6)
}

# Counts of things annotated with particular tissues
# (including things not significant.)
tissue.names = sort(unique(ao$group.name))
tissue.lineage = matrix(0, nrow=23, ncol=length(tissue.names))
rownames(tissue.lineage) = cells.23
colnames(tissue.lineage) = tissue.names
for(cell in cells.23) {
  g = lin.enriched[ lin.enriched$set == paste0(cell, "_enriched") , "gene" ]
  at = ao[ ao$gene %in% g, "group.name" ]

  for(tissue in tissue.names) {
    tissue.lineage[cell,tissue] = sum(at==tissue)
  }
}

interior.cells =
  names(which(apply(cell.lineage.matrix[,cells.23] > 0, 1, any)))

pdf("git/sort_paper/plot/lineage/lineageEnriched.pdf",
  width=15, height=7)
r = data.frame(cell = as.character(rownames(cell.time.on.off)),
  time.1=cell.time.on.off$on, time.2=cell.time.on.off$off,
  stringsAsFactors=FALSE)
r = r[ r$cell %in% interior.cells , ]
r$col = "#c0c0c0"
r$x = number.lineage.columns(lin)[ r$cell ]
rownames(r) = r$cell

# plot the tree
plot.segments(r, "Anatomy annotations in lineages",
  root="P0", c(0,133), lwd=3, yaxt="n",
  int.n.to.label = cells.23, srt=90)

# rectangle to mark known annotation
rect(-10000, 123, 10000, 150, lwd=0, col="lightgrey")

pie.colors =
  hsv(0:(length(tissue.names)-1) / length(tissue.names), 0.7, 0.9)

for(cell in cells.23) {
  tl = tissue.lineage[ cell , ]

  # first, the measured counts
  if (sum(tl) > 0) {
    floating.pie(r[cell,"x"], r[cell,"time.2"],
      tl + 1e-6, col = pie.colors, radius = 9)
  }

  # then, the known tissue types per lineage (from annotation)
  floating.pie(r[cell,"x"], 130,
    known.tissue.count(cell) + 1e-6,
    col=pie.colors, radius = 9)
}

axis(2, at=c(100, 130), labels=c("predicted", "known"), tick = FALSE)

legend("topleft", legend=tissue.names, fill = pie.colors, cex=0.84)

dev.off()

