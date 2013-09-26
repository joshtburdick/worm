# Computes upstream regions.
# For now, using the WS220 build from Cufflinks.
# Deprecated -- switching to using bedtools.
# Nope, this may actually be used.

out.file = "git/tf/motif/motifCount/upstream_liftOver_WS220.bed"

# sizes of chromosomes (used to prevent regions that are off the end)
chr.sizes = {
  s = read.table("/home/jburdick/data/seq/Caenorhabditis_elegans/Ensembl/WS220/Sequence/WholeGenomeFasta/genome.fa.fai")
  r = s[,2]
  names(r) = s[,1]
  r
}

# base annotation table
g = read.table("git/unmix/seq/quant/geneBounds_liftOver_WS220.bed",
  sep="\t", as.is=TRUE)
colnames(g) = c("chr", "a", "b", "gene.id", "score", "strand")

# Computes upstream intergenic regions.
# Args:
#   g - a BED file, as a data.frame
#   min.dist - the minimum upstream distance to include
#   max.dist - the maximum upstream distance to include
# Returns: g, with extra columns "a1" and "b1", giving
#   the coordinates of the relevant upstream region.
compute.upstream.intergenic = function(g, min.dist, max.dist) {

  g$a1 = NA
  g$b1 = NA
cat("\n")
  for(i in 1:nrow(g)) {
cat("\b\b\b\b\b\b\b", i)
    # XXX this is a pretty inefficient way of doing this.
    if (g$strand[i] == "+") {
      g$a1[i] =
        max(c(1, g$b[ g$chr==g$chr[i] & g$b < g$a[i] ]))      
    }

    if (g$strand[i] == "-") {
      g$b1[i] =
        min(c(chr.sizes[g$chr[i] ] - 1, g$a[ g$chr==g$chr[i] & g$a > g$b[i] ]))      
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


  stopifnot(all(g$a1 <= g$b1))

  g
}

upstream.region = compute.upstream.intergenic(g, 500, 5000)
upstream.region$a = upstream.region$a1
upstream.region$b = upstream.region$b1
upstream.region = upstream.region[,c(1:6)]

if (TRUE) {
write.table(upstream.region,
  file=out.file,
  sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

system(paste("bedSort", out.file, out.file))
}

