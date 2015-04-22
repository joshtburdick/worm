# Looks for unannotated transcription in Cufflinks' de novo output.

source("git/utils.r")

cuff.result.dir = "/media/jburdick/disk2/jburdick/cufflinks_denovo/"

output.dir = "git/sort_paper/FACS/otherTranscription/coverageAndLength/"

system(paste("mkdir -p", output.dir))

# Gets statistics about transcripts (annotated or not) from one experiment.
# Args:
#   cuff.dir - directory containing Cufflinks' output
# Returns: data frame of statistics about transcripts
transcript.stats = function(cuff.dir) {
  genes = read.table(paste0(cuff.result.dir, "/", cuff.dir, "/genes.fpkm_tracking"),
    header=TRUE, sep="\t", as.is=TRUE)
  isoforms = read.table(paste0(cuff.result.dir, "/", cuff.dir, "/isoforms.fpkm_tracking"),
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
# Side effects: writes a GTF file
write.gtf.subset = function(genes, in.file, out.file) {
  a = read.tsv(in.file)

  # parse out the gene ID
  p = 'gene_id "([^"]+)";'
  a$gene.id = sapply(regmatches(x1, regexec(p, x1)), function(a) a[2])

  a
#  write.tsv(r, out.file)
}

# a = transcript.stats("Murray_050912/cehm36p")





# Plots the coverage and length of transcripts, from Cufflinks or otherwise.
plot.transcript.stats = function(r, name) {
  png(paste0(output.dir, "/", name, ".png"), width=1000, height=1000)
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

plot.transcript.stats(a, "ceh-36 (+)")

