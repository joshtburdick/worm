# A sort matrix, corresponding to the groups of cells
# measured by Spencer et al.

library(Matrix)

source("git/utils.r")

load("R/lineage/tree_utils.Rdata")
load("R/lineage/cell.names.Rdata")
load("git/data/wormbase/anatomy.ontology.group.Rdata")

ao.name.to.term = c(tapply(names(ao.term.to.name), ao.term.to.name, function(a) a[1]))

# tissue.by.lineage = read.table("data/worm/TissuesPerCell.tsv",
#   sep="\t", header=TRUE, quote="", row.names=1, as.is=TRUE)
# nervous system cells, obtained a different way
# nerve.cells = rownames(tissue.by.lineage)[
#   tissue.by.lineage[,"Tissue"] %in% c("Glia","Nervous") ]

# Given the name of something, gets the associated cell names.
# Args:
#   ao.name - the name of an anatomy ontology
# Returns: a boolean vector, indexed by lineage terms
anatomy.term.to.lineage = function(ao.name) {

  # extract cell linage names
  ao.term = ao.name.to.term[ ao.name ]
  stopifnot(!is.na(ao.term))

  # get all contained anatomy terms
  ao.term = intersect(ao.term, colnames(ao.subterm))
  ao.term.1 = ao.term.to.name[
    colnames(ao.subterm)[ ao.subterm[ , ao.term ] ] ]

  # get just names of nuclei
#  ao.term.1 = grep(" nucleus$", ao.term.1, value=TRUE) 
  ao.term.1 = sub(" nucleus$", "", ao.term.1)
  ao.term.1 = sub("\\.", "", ao.term.1)

  ao.cell = lin.node.names %in% ao.term.1
  names(ao.cell) = lin.node.names

  ao.cell
}

# names of lineage terms to include (this is only the
# embryonic cell lines)
spencer.ao.terms = c("BAG",
  "germline precursor cell",
  "AVA",
  "GABAergic neuron",
  "body wall muscle cell",     # or "body wall musculature"
  "coelomocyte",
  "dopaminergic neuron",
  "intestinal cell",     # or "intestine"  
  "nervous system",
  "hypodermis",
  "AVE",
  "pharyngeal muscle cell")

spencer.m = t(sapply(spencer.ao.terms, anatomy.term.to.lineage))

# A-class motor neurons (unc-4 also is on elsewhere, esp. later)
spencer.m = rbind(spencer.m, "A-class motor.neuron" =
  anatomy.term.to.lineage("DA neuron") | anatomy.term.to.lineage("VA neuron"))

# Adds in cells, both of whose daughter cells are on.
# Args:
#   x - a boolean vector, indexed by cell name
# Returns: x, with parent lineages filled in.
fill.in.lineage = function(x) {
  x.names = names(x)[x]
  p1 = lin.triples[ lin.triples$daughter.1 %in% x.names , "parent" ]
  p2 = lin.triples[ lin.triples$daughter.2 %in% x.names , "parent" ]
  p = intersect(p1, p2)
  x1 = lin.node.names %in% union(x.names, p)
  names(x1) = lin.node.names
  if (all(x == x1))
    x
  else
    fill.in.lineage(x1)
}

# do this for all these sets of cells
spencer.m = t(apply(spencer.m, 1, fill.in.lineage))

spencer.m = spencer.m[ order(apply(spencer.m, 1, sum)) , ]

write.tsv(spencer.m, "git/sort_paper/validation/spencerSortMatrix.tsv")

# pdf("git/sort_paper/validation/spencerSortMatrix.pdf", width=10, height=5)
png("git/sort_paper/validation/spencerSortMatrix.png", width=1000, height=500)

par(mar=c(0.1,11,1,0.1)+0.1)

image(t(spencer.m),
  useRaster=TRUE,
  xaxt="n", yaxt="n", zlim=c(0,1),
  col=hsv(1/6, 1, 0:128/128))

axis(2, at=c(0:12)/12,
  labels=rownames(spencer.m),
  tick=FALSE, las=2)

dev.off()

