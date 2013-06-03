# Makes a BED file using "common" gene names.
# This only includes the overall boundaries,
# ignoring exon structure.

# Gets transcript bounds from a GFF file.
get.gff.transcript.bounds = function(gff.file) {
  a = read.table(gff.file, sep="\t", stringsAsFactors=FALSE)
  a = a[ a[,3] == "transcript" , ]
  n = a[,9]
  r = data.frame(chr = a[,1], start = a[,4], end = a[,5],
    name = sub(";$", "", sub(".*transcript_id ", "", a[,9])),
    score = 0, strand = a[,7],
    stringsAsFactors = FALSE)
  rownames(r) = r$name
  r
}

# Forwarding declaration.
get.lincrna.bounds = get.gff.transcript.bounds

ncrnas = rbind(
  get.lincrna.bounds("~/gcb/data/expression/lincRNA/ancRNAs_W3PSeq3_ce6.gtf"),
    get.lincrna.bounds("~/gcb/data/expression/lincRNA/lincRNAs_W3PSeq3_ce6.gtf"))

# Get gene bounds from (I think this is .psl format?)
get.transcript.bounds = function(f) {
  g = read.table(f, sep="\t", as.is=TRUE)

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

#  b = rbind(b, c("chrI", 15058993, 15072423, "ribosomal_RNA", "0", "+"))
  b
}

# Currently using ce6.
gene.bounds = get.transcript.bounds(gzfile("~/data/seq/Caenorhabditis_elegans/UCSC/ce6/Annotation/Genes/refFlat.txt.gz"))

b = rbind(gene.bounds, ncrnas)

# tack on "ribosomal RNAs"
b = rbind(b, c("chrI", 15058993, 15072423, "ribosomal_RNA", "0", "+"))

# Gets bounds for WS220 (ce10) genes.
get.gene.bounds.ws220 = function() {

  gene.bounds = get.transcript.bounds(
    "/home/jburdick/data/seq/Caenorhabditis_elegans/Ensembl/WS220/Annotation/Genes/refFlat.txt.gz")
  gene.bounds = gene.bounds[ -grep("^21ur-", gene.bounds$name) , ]
  ncrna.bounds = rbind(get.gff.transcript.bounds(
    "/home/jburdick/gcb/data/expression/lincRNA/ancRNAs_W3PSeq3_ce10.gtf"),
    get.gff.transcript.bounds(
    "/home/jburdick/gcb/data/expression/lincRNA/lincRNAs_W3PSeq3_ce10.gtf"))
  ncrna.bounds$chr = sub("chr", "", ncrna.bounds$chr)
  rbind(gene.bounds, ncrna.bounds,
    c("I", 15059200, 15066508, "ribosomal_RNA", "0", "+"))
}

if (TRUE) {
  write.table(b, file="git/unmix/seq/quant/geneBounds.tsv",
    sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

  # write out gene bounds on the opposite strand (which, confusingly,
  # are the bounds used for "sense" gene quantification)
  b.flip = b
  b.flip$strand = ifelse(b$strand == "-", "+", "-")
  write.table(b.flip, file="git/unmix/seq/quant/geneBounds_flip.tsv",
    sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

  # similarly for WS220
  b.ws220 = get.gene.bounds.ws220()

  write.table(b.ws220,
    file="git/unmix/seq/quant/geneBounds_WS220.tsv",
    sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  b.ws220.flip = b.ws220
  b.ws220.flip$strand = ifelse(b.ws220.flip$strand == "-", "+", "-")
  write.table(b.ws220.flip,
    file="git/unmix/seq/quant/geneBounds_WS220_flip.tsv",
    sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

a = get.gene.bounds.ws220()


