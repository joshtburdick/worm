# Computes upstream regions.
# For now, using the WS220 build from Cufflinks.
# Deprecated -- switching to using bedtools.
# Nope, this may actually be used

out.file = "git/tf/motif/motifCount/upstreamRegionsWS220_5kb_nogenes.bed"

upstream.dist = 5000

# base annotation table
g = read.table("git/unmix/seq/quant/geneBounds_WS220.bed",
  sep="\t", as.is=TRUE)
colnames(g) = c("chr", "a", "b", "gene.id", "score", "strand")


# Computes upstream intergenic regions.
# Args:
#   g - a BED file, as a data.frame
#   max.dist - the maximum upstream distance to include
# Returns: g, with extra columns "a1" and "b1", giving
#   the coordinates of the relevant upstream region.
compute.upstream.intergenic = function(g, max.dist) {

  g$a1 = NA
  g$b1 = NA
cat("\n")
  for(i in 1:nrow(g)) {
cat("\b\b\b\b\b\b\b", i)
    # XXX this is a pretty inefficient way of doing this.
    if (g$strand[i] == "+") {
      g$a1[i] =
        max(c(0, g$b[ g$chr==g$chr[i] & g$b < g$a[i] ]))      
    }

    if (g$strand[i] == "-") {
      g$b1[i] =
        min(c(1e10, g$a[ g$chr==g$chr[i] & g$a > g$b[i] ]))      
    }

  }
cat("\n")
  # tack on other end
  g$a1[ g$strand == "-" ] = g$b[ g$strand == "-" ]
  g$b1[ g$strand == "+" ] = g$a[ g$strand == "+" ]

  # trim length to no more than max.dist
  i = g$strand=="+" & (g$b1 - g$a1) > max.dist
  g$a1[ i ] = g$a[ i ] - max.dist
  i = g$strand=="-" & (g$b1 - g$a1) > max.dist
  g$b1[ i ] = g$b[ i ] + max.dist

  # force this to be positive
  g$a1[ g$a1 < 0 ] = 0

  g
}

# compute upstream regions
upstream.region = data.frame(chr=g$chr,
  a=ifelse(g$strand=="+", g$a - upstream.dist, g$b),
  b=ifelse(g$strand=="+", g$a, g$b + upstream.dist),
  gene=g$g,
  score=0,
  strand=g$strand)
upstream.region$a[ upstream.region$a < 1 ] = 1

upstream.region = compute.upstream.intergenic(g, 5000)
upstream.region$a = upstream.region$a1
upstream.region$b = upstream.region$b1
upstream.region = upstream.region[,c(1:6)]

write.table(upstream.region,
  file=out.file,
  sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

system(paste("bedSort", out.file, out.file))

