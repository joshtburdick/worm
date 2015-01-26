# Plots coverage, as a function of distance (in transcript space)
# from the 3-prime end.
# Other possible ways to do this:
# - SeqMonk
# - RseQC

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

# Gets the 3'-most end of a vector.
# Args:
#   dist.3p - the length of 3' portion to include
#   x - an Rle object (or vector), of the dist.3p"
# Returns: the "dist.3p"-most entries of x (left-padded
#   with 0's if necessary)
get.3p.end = function(dist.3p, x) {
  n = length(x)

  # is this long enough?
  if (n >= dist.3p)
    # yes: return the end of it
    x[ (n+1-dist.3p) : n ]
  else
    # no: left-pad it with 0's
    c( Rle(0, dist.3p - n), x)
}


# Gets the 3' coverage

# Args:


# Returns: list containing
#   coverage.3p - vector of coverage at the 3' ends
#   gene.length - vector of the distribution of the lengths
#     of the exonic regions
get.3p.coverage = function(r, exon.loc, gene.loc, dist.3p) {

  # accumulators for the coverage stats, and lengths of genes
  coverage.3p = rep(0, dist.3p)
  gene.length = c()

  # loop through each strand and chromosome
  for(strand in c("+", "-"))
    for(chr in setdiff(names(ov.plus), "MtDNA")) {
      write.status(c(chr, strand))
      
      exon.loc.1 = exon.loc[
        exon.loc$space==chr & exon.loc$strand==strand , ]
      gene.loc.1 = gene.loc[
        gene.loc$space==chr & gene.loc$strand==strand , ]
      r1 = r[[ chr ]]

      # We need to assign each exonic region to a gene.
      # On the "+" strand, we want to pick the first gene
      # overlapping a given exon, since presumably it has the
      # nearest 3' end (and vice versa on the "-" strand).
      exon.gene.overlap = findOverlaps(exon.loc.1, gene.loc.1,
        select = ifelse(strand=="+", "first", "last"))

      # get indices of exonic region, grouped by gene
      i = by(1:length(exon.gene.overlap), exon.gene.overlap, c)

      e1 = exon.loc.1$ranges

      # loop through the grouped exonic regions
      for(j in i[1:10]) {

        exon.regions = e1[j,]

        c1 = r1[ exon.regions ]
        if (strand=="-") c1 = rev(c1)
        coverage.3p = coverage.3p + get.3p.end(c1, dist.3p)

        gene.length = c(gene.length, sum(width(exon.regions)))
      }

  }

  list(coverage.3p = coverage.3p, gene.length = gene.length)
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

foo = get.3p.coverage(r, exon.loc, gene.loc, 10000)


