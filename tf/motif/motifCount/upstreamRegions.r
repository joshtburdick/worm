# Computes upstream regions.
# For now, using the WS220 build from Cufflinks.

source("git/data/name_convert.r")

out.path = "git/tf/motif/motifCount/regions/"

system(paste("mkdir -p", out.path))

# classes of genes (coding, ncRNA, etc.)
gene.class = read.table(gzfile(
  "git/data/seq/gene_class.tsv.gz"), as.is=TRUE)
colnames(gene.class) = c("gene", "class")
rownames(gene.class) = gene.class$gene
gene.class = rename.gene.names(gene.class)

# sizes of chromosomes (used to prevent regions that are off the end)
chr.sizes = {
  s = read.table("/home/jburdick/data/seq/Caenorhabditis_elegans/Ensembl/WS220/Sequence/WholeGenomeFasta/genome.fa.fai")
  r = s[,2]
  names(r) = s[,1]
  # XXX
  names(r)[ names(r) == "MT" ] = "MtDNA"
  r
}

# base annotation table
g = read.table("git/data/seq/merged_genes_WS220.bed",
  sep="\t", as.is=TRUE)
colnames(g) = c("chr", "a", "b", "gene.id", "score", "strand")

# Computes upstream intergenic regions.
# Args:
#   g - a BED file, as a data.frame
#   g1 - a similar BED file, of genes to use to define upstream regions
#   min.dist - the minimum upstream distance to include
#   max.dist - the maximum upstream distance to include
# Returns: g, with extra columns "a1" and "b1", giving
#   the coordinates of the relevant upstream region.
compute.upstream.intergenic = function(g, g1, min.dist, max.dist) {

  g$a1 = NA
  g$b1 = NA
cat("\n")
  for(i in 1:nrow(g)) {
cat("\b\b\b\b\b\b\b", i)
    # XXX this is a pretty inefficient way of doing this
    if (g$strand[i] == "+") {
      # if any gene overlaps gene start, set its upstream size to 0
      if (any((g1$chr==g$chr[i]) & (g1$a < g$a[i]) & (g1$b > g$a[i]))) {
        g$a1[i] = g$a[i]
      }
      # otherwise, find distance to nearest gene
      else {
        g$a1[i] = max(c(1,
          g1$b[ g1$chr==g$chr[i] & g1$b < g$a[i] ]))      
      }
    }

    if (g$strand[i] == "-") {
      # symmetric case
      if (any((g1$chr==g$chr[i]) & (g1$a < g$b[i]) & (g1$b > g$b[i]))) {
        g$b1[i] = g$b[i]
      }
      else {
        g$b1[i] = min(c(chr.sizes[g$chr[i] ] - 1,
          g1$a[ g1$chr==g$chr[i] & g1$a > g$b[i] ]))      
      }
    }
  }
cat("\n")

  # tack on other end
  g$a1[ g$strand == "-" ] = g$b[ g$strand == "-" ]
  g$b1[ g$strand == "+" ] = g$a[ g$strand == "+" ]

  # trim length (XXX slightly ugly)
  # min
  i = g$strand=="+" & (g$b1 - g$a1) < min.dist
  g$a1[ i ] = g$a[ i ] - min.dist
  i = g$strand=="-" & (g$b1 - g$a1) < min.dist
  g$b1[ i ] = g$b[ i ] + min.dist

  # max
  i = g$strand=="+" & (g$b1 - g$a1) > max.dist
  g$a1[ i ] = g$a[ i ] - max.dist
  i = g$strand=="-" & (g$b1 - g$a1) > max.dist
  g$b1[ i ] = g$b[ i ] + max.dist

  # trim these so they are still on the chromosomes
  g$a1[ g$a1 < 0 ] = 0
  g$b1[ g$b1 <= g$a1 ] = g$a1[ g$b1 <= g$a1 ] + 1
  s = chr.sizes[ g$chr ] - 1
  g$b1[ g$b1 >= s ] = s[ g$b1 >= s ]
  g$a1[ g$a1 >= g$b1 ] = g$b1[ g$a1 >= g$b1 ] - 1

#  stopifnot(all(g$a1 <= g$b1))

  g$a = g$a1
  g$b = g$b1
  g[ , c(1:6) ]
}

# coding genes (used to mark boundaries of regions
# upstream of genes)
coding.genes = g[ g$gene.id %in%
  rownames(gene.class)[gene.class$class == "protein_coding" ] , ]



for(upstream.size in c(1:5)) {
  cat("\n", upstream.size, "\n")
  upstream.region = compute.upstream.intergenic(
    g, coding.genes, 500, 1000 * upstream.size)
  out.file = paste(out.path, "/upstream_", upstream.size,
    "kb_WS220.bed", sep="")
  write.table(upstream.region,
    file=out.file,
    sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  system(paste("bedSort", out.file, out.file))
}

