# Looks for unannotated transcription in Cufflinks' de novo output.

source("git/utils.r")
source("git/data/gtf.r")

cuff.result.dir = "/media/jburdick/disk2/jburdick/cufflinks_denovo/"

output.dir = "git/sort_paper/FACS/otherTranscription/"

system(paste0("mkdir -p ", output.dir, "coverageAndLength"))

# Gets statistics about transcripts (annotated or not) from one experiment.
# Args:
#   cuff.dir - directory containing Cufflinks' output
# Returns: data frame of statistics about transcripts
transcript.stats = function(cuff.dir) {
  genes = read.table(paste0(cuff.dir, "/genes.fpkm_tracking"),
    header=TRUE, sep="\t", as.is=TRUE)
  isoforms = read.table(paste0(cuff.dir, "/isoforms.fpkm_tracking"),
    header=TRUE, sep="\t", as.is=TRUE)
  gene.length = tapply(isoforms$length, isoforms$gene_id, max)
  genes$length = gene.length[ genes$gene_id ]
  genes$cuff = grepl("CUFF.", genes$gene_id)
  genes
}

# Writes out a subset of a GTF file, containing just transcripts
# called as "novel."
# Args:
#   genes - list of genes to include
#   in.file, out.file - GTF file to read and write, respectively
# Returns: the parsed file of results.
write.gtf.subset = function(genes, in.file, name) {

  # read in just the new transcripts which didn't overlap
  in.file = paste0(cuff.result.dir, "Murray_050912/cehm36p/transcripts.gtf")
  a = read.table(sep="\t", quote="", as.is=TRUE,
    pipe(paste0("bedtools subtract -s -A ",
    " -a ", f,
    " -b git/data/seq/merged_genes_WS220.bed ")))

  # parse out the gene ID
  p = 'gene_id "([^"]+)";'
  x1 = a[,9]
  a$gene.id = sapply(regmatches(x1, regexec(p, x1)), function(m) m[2])

  # limit to the subset
  a = a[ a$gene.id %in% genes , ]

#  write.table(a[ , c(1:9) ], paste0(output.dir, 
}



# Plots the coverage and length of transcripts, from Cufflinks or otherwise.
# Args:
#   r - transcript stats
#   name - name to use for output file
plot.transcript.stats = function(r, name) {
  png(paste0(output.dir, "/coverageAndLength/", name, ".png"), width=1000, height=1000)
  xlim=c(0, max(log10(r$length)))
  ylim=c(0, max(log2(1+r$FPKM)))
  plot(log10(r$length), log2(1+r$FPKM), xlim=xlim, ylim=ylim,
    main=name, xlab="log10(length)", ylab="log2(1+FPKM)",
    pch=20, cex=0.7, col="#00000040")
  par(new=TRUE)
  r1 = r[r$cuff,]
  plot(log10(r1$length), log2(1+r1$FPKM), xlim=xlim, ylim=ylim,
    main="", xlab="", ylab="",
    pch=20, cex=0.7, col="#000000a0")

  dev.off()
}

# Summarizes transcripts from one experiment.
# Args:
#   input.dir - directory containing Cufflinks' output
#   name - name to use for output files
# Side effects: writes out files.
summarize.file = function(input.dir, name) {
  r = transcript.stats(input.dir)

  # read in just the new transcripts which didn't overlap
  gtf = read.table(sep="\t", quote="", as.is=TRUE,
    pipe(paste0("bedtools subtract -s -A ",
    " -a ", paste0(input.dir, "/transcripts.gtf"),
    " -b git/data/seq/merged_genes_WS220.bed ")))



  plot.transcript.stats(r, name)
}

summarize.file(paste0(cuff.result.dir, "Murray_050912/cehm36p"), "ceh-36 (+)")

