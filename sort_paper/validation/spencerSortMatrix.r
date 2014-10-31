# A sort matrix, corresponding to the groups of cells
# measured by Spencer et al.

library(Matrix)

load("R/lineage/tree_utils.Rdata")
load("R/lineage/cell.names.Rdata")
load("git/data/wormbase/anatomy.ontology.group.Rdata")

ao.name.to.term = c(tapply(names(ao.term.to.name), ao.term.to.name, function(a) a[1]))

# tissue.by.lineage = read.table("data/worm/TissuesPerCell.tsv",
#   sep="\t", header=TRUE, quote="", row.names=1, as.is=TRUE)

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
  "body wall muscle cell",       # or "body wall musculature"
  "coelomocyte",
  "dopaminergic neuron",
  "intestinal cell",       # or "intestine"  
  "nervous system",
  "hypodermis",
  "AVE",
  "pharyngeal muscle cell")

# A-class motor neurons (unc-4 also is on elsewhere, esp. later)
# "DA neuron", "VA neuron",

# nervous system cells, obtained a different way
nerve.cells = rownames(tissue.by.lineage)[
  tissue.by.lineage[,"Tissue"] %in% c("Glia","Nervous") ]






