# Computes distance from the 3' end for alternative splicing
# events.
# XXX not very fast

source("git/utils.r")

# gene boundaries
gene.bounds = read.table("git/data/seq/merged_genes_WS220.bed", as.is=TRUE)
colnames(gene.bounds) = c("chr", "a", "b", "name", "score", "strand")
rownames(gene.bounds) = gene.bounds$name

# all exons of each gene, merged together
exons = read.table("git/data/seq/merged_exons_WS220.bed", as.is=TRUE)
colnames(exons) = c("chr", "a", "b", "name", "score", "strand")
exons$size = exons$b - exons$a

# Finds distance from the 3' end of a gene to some
# location (e.g. a splicing event.) This distance
# includes all annotated exons.
# Args: r - a data frame with columns including
#     gene - gene name
#     chrom, pos, strand - genomic location
# Returns: the same data frame, but with a new column
#   dist.to.3p - the total exonic distance from the 3'
#     end to the genomic location (not including any
#     overlapping exons)
get.distance.to.3.prime = function(r) {
  r$dist.to.3p = NA

  for(i in 1:nrow(r)) {
    write.status(paste(i, "/", nrow(r)))

    if (r[i,"gene"] %in% rownames(gene.bounds)) {
      g = gene.bounds[ r[i,"gene"] , ]
      
      if (g$strand == "+") {
        bounds = c(r[i,"pos"], g$b)
      }
      else if (r[i,"strand"] == "-") {
        bounds = c(g$a, r[i,"pos"])
      }

      e1 = exons[ exons$chr == g$chr &
        exons$a >= bounds[1] &
        exons$b <= bounds[2] &
        exons$strand == g$strand , ]
      r[i,"dist.to.3p"] = sum(e1$size)
    }
  }
  cat("\n")

  r
}

event.names = read.table("git/unmix/seq/quant/MISO/geneAndAS.tsv",
  header=TRUE, as.is=TRUE)
event.names$pos = ifelse(event.names$strand=="+", event.names$end, event.names$start)
event.names = get.distance.to.3.prime(event.names)

event.names = event.names[ , names(event.names) != "pos" ]

write.tsv(event.names,
  "git/sort_paper/splicing/MISO/asEventAnnotation.tsv")

