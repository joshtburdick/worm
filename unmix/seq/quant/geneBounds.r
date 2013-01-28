# Makes a BED file using "common" gene names.
# This only includes the overall boundaries,
# ignoring exon structure.

# Gets bounds of lincRNAs, from Bartel lab GFF files.
get.lincrna.bounds = function(gff.file) {
  a = read.table(gff.file, sep="\t")
  a = a[ a[,3] == "transcript" , ]
  n = a[,9]
  r = data.frame(chr = a[,1], start = a[,4], end = a[,5],
    name = sub(";$", "", sub(".*transcript_id ", "", a[,9])),
    score = 0, strand = a[,7], stringsAsFactors = FALSE)
  rownames(r) = r$name
  r
}

ncrnas = rbind(
  get.lincrna.bounds("~/gcb/data/expression/lincRNA/ancRNAs_W3PSeq3_ce6.gtf"),
    get.lincrna.bounds("~/gcb/data/expression/lincRNA/lincRNAs_W3PSeq3_ce6.gtf"))

# Get the UCSC annotations.
get.annotations = function() {
  g = read.table(gzfile("~/data/seq/Caenorhabditis_elegans/UCSC/ce6/Annotation/Genes/refFlat.txt.gz"),
    sep="\t", as.is=TRUE)

  colnames(g) = c("name", "refseq.id", "chr", "strand", "tx.start",
    "tx.end", "cds.start", "cds.end",
    "exon.start", "exon.end")
  g = g[ g$name != "", ]

  # utility checking that something is always the same
  all.same = function(x)
    { x1 = as.vector(x); if (all(x1==x1[1])) x1[1] else NA }

  g.chr = c(by(g$chr, g$name, all.same))
  g.strand = c(by(g$strand, g$name, all.same))
  g.start = c(by(g$tx.start, g$name, min))
  g.end = c(by(g$tx.end, g$name, max))

  names = unique(g$name)

  b = data.frame(chr = g.chr[names],start=g.start[names], end=g.end[names],
    name = names, score = 0, strand = g.strand[names], stringsAsFactors=FALSE)

  b = b[ !is.na(b$chr) , ]
  b = b[ !is.na(b$strand) , ]

  # tack on "ribosomal RNAs"
  b = rbind(b, c("chrI", 15058993, 15072423, "ribosomal_RNA", "0", "+"))
  b
}

gene.bounds = get.annotations()
b = rbind(gene.bounds, ncrnas)

if (TRUE) {
write.table(b, file="git/unmix/seq/quant/geneBounds.tsv",
  sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

  # write out gene bounds on the opposite strand (which, confusingly,
  # are the bounds used for "sense" gene quantification)
  b.flip = b
  b.flip$strand = ifelse(b$strand == "-", "+", "-")
  write.table(b.flip, file="git/unmix/seq/quant/geneBounds_flip.tsv",
    sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}


