# Computes upstream regions.
# For now, using the WS220 build from Cufflinks.

out.file = "git/tf/motif/motifCount/upstreamRegionsWS220_10kb.bed"

upstream.dist = 10000

# base annotation table
g1 = read.table("~/data/seq/Caenorhabditis_elegans/Ensembl/WS220/Annotation/Genes/refGene.txt",
  as.is=TRUE)
g = data.frame(gene.id=g1[,2], gene=g1[,13],
  chr=g1[,3], a=g1[,5], b=g1[,6], strand=g1[,4], stringsAsFactors=FALSE)
g = g[ substr(g$gene, 1, 5) != "21ur-" , ]

# for now, we only want to include genes in this list
geneBounds = read.table("git/unmix/seq/quant/geneBounds.tsv", as.is=TRUE)
geneList = geneBounds[,4]

g$g = NA
g$g[ g$gene %in% geneList ] = g$gene[ g$gene %in% geneList ]
g$g[ g$gene.id %in% geneList ] = g$gene.id[ g$gene.id %in% geneList ]
g = g[ !is.na(g$g) , ]
g = g[ !duplicated(g$g) , ]

# compute upstream regions
upstream.region = data.frame(chr=g$chr,
  a=ifelse(g$strand=="+", g$a - upstream.dist, g$b),
  b=ifelse(g$strand=="+", g$a, g$b + upstream.dist),
  gene=g$g,
  score=0,
  strand=g$strand)
upstream.region$a[ upstream.region$a < 1 ] = 1

write.table(upstream.region,
  file=out.file,
  sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

system(paste("bedSort", out.file, out.file))

