# Plots coverage, as a function of distance (in transcript space)
# from the 3-prime end.

library("rtracklayer")

source("git/utils.r")

bw.file = "/murrlab/seq/igv/expression/embryo_FACS//wig/irx-1 (+).bw"

# FIXME: ideally the set of genes should be non-overlapping
# (as some of these overlap)
gene.loc = import("git/data/seq/merged_genes_WS220.bed")

# note that this is "merged exonic regions", not "exons" per se
exon.loc = import("git/data/seq/merged_exons_WS220.bed")

# a gene (possibly) for testing
test.gene = list(g1 = exon.loc["I"][7:11,])


# We need to assign each exonic region to a gene.
# On the "+" strand, we want to pick the first gene
# overlapping a given exon, since presumably it has the
# nearest 3' end (and vice versa on the "-" strand).
exon.plus = exon.loc[exon.loc$strand=="+",]
exon.minus = exon.loc[exon.loc$strand=="-",]
ov.plus = findOverlaps(exon.plus,
  gene.loc[ gene.loc$strand == "+" , ], select="first")
ov.minus = findOverlaps(exon.minus,
  gene.loc[ gene.loc$strand == "-" , ], select="last")

x = NULL
for(chr in setdiff(names(ov.plus), "MtDNA")) {
#  chr = "I"
  write.status(chr)

  # get indices of exonic region, grouped by gene
  i = by(1:length(ov.plus[[chr]]), ov.plus[[chr]], c)

  e1 = exon.plus[chr]$ranges
  for(j in sample(i, 10))
    x = c(x, e1[j,])

}

# by(1:length(ov.minus[[chr]]), ov.minus[[chr]], c)


# Construct a list of all the exonic regions for each gene.
# Args:
#   gene.loc - the overall gene bounds
#   exon.loc - the non-overlapping exon regions
# Returns: list, indexed by gene name, of SimpleRleList
#   objects.
get.regions.per.gene = function(gene.loc, exon.loc) {







}




# Finds total coverage in all the exons of a gene.
# Possibly deprecated.
# Args:
#   r - a SimpleRleList object
#   exon.regions - a RangedData object, representing the exonic
#     regions of one gene (not the exons, but rather
#     those regions merged)
# Returns: an Rle object of coverage at each exonic region
coverage.relative.to.3p = function(r, exon.regions) {
  g = exon.regions[[1]][1]

  stopifnot(all(g$space == g$space[1]))
  stopifnot(all(g$strand == g$strand[1]))

  # the space (in this case, chromosome) to get data from
  sp = as.character(g$space[1])

  # get concatenated coverage for all these regions
  # this is smoother than I'd expected!
  a = r[[sp]][ g$ranges ]

  # reverse coverage if strand is "-"
  if (g$strand[1] == "-") a = rev(a)

  a
}


# gets a SimpleRleList, with coverage for each chromosome
r = coverage(import.bw(bw.file))



# foo = coverage.relative.to.3p(r, test.gene)



