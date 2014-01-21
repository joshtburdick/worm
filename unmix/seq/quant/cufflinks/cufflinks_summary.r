# Summarizes results from Cufflinks analysis.
# Something of an oversimplification; doesn't include

source("git/utils.r")

input.path = "git/unmix/seq/quant/cufflinks/cufflinks_out_20140121/"
output.path = "git/unmix/seq/quant/cufflinks/cufflinks_summary_20140121/"

# Gets results from one directory.
get.cufflinks.results = function(f) {

  genes = read.table(
    gzfile(paste(f, "/genes.fpkm_tracking.gz", sep="")),
    sep="\t", header=TRUE, as.is=TRUE)
  genes = genes[ ! grepl("21ur-", genes$gene_short_name) , ]
  genes = genes[ , c(1,5,10,13) ]
  genes[ genes$FPKM_status != "OK" , "FPKM" ] = 0

  # XXX hack to default to using a "human readable" name
  genes$id = ifelse(genes$gene_short_name == "-",
    genes$tracking_id, genes$gene_short_name)

  # XXX this is to handle the very few cases in which a
  # gene name is repeated
  gene.fpkm = c(by(genes$FPKM, genes$id, sum))

  isoforms = read.table(
    gzfile(paste(f, "/isoforms.fpkm_tracking.gz", sep="")),
    sep="\t", header=TRUE, as.is=TRUE)
  isoforms = isoforms[ ! grepl("21ur-", genes$gene_short_name) , ]
  isoforms = isoforms[ , c(1,5,10,13) ]

  isoform.fpkm = c(by(isoforms$FPKM, isoforms$tracking_id, sum))

  list(genes = genes, isoforms = isoforms,
    gene.fpkm = gene.fpkm, isoform.fpkm = isoform.fpkm)
}


# Converts all results in a directory of Cufflinks results.
convert.dir = function(input.dir, output.basename) {
  f = list.files(input.dir)

  a1 = get.cufflinks.results(paste(input.dir, "/", f[1], sep=""))

  gene.fpkm = matrix(NA, nrow=length(a1$gene.fpkm), ncol=length(f))
  dimnames(gene.fpkm) = list(sort(names(a1$gene.fpkm)), f)
  isoform.fpkm = matrix(NA, nrow=length(a1$isoform.fpkm), ncol=length(f))
  dimnames(isoform.fpkm) = list(sort(names(a1$isoform.fpkm)), f)

  for(j in 1:length(f)) {
    cat(j, backspace.string, f[j])
    r = get.cufflinks.results(paste(input.dir, "/", f[j], sep=""))
    gene.fpkm[,j] = r$gene.fpkm[ rownames(gene.fpkm) ]
    isoform.fpkm[,j] = r$isoform.fpkm[ rownames(isoform.fpkm) ]
  }

  write.tsv(gene.fpkm,
    file=gzfile(paste(output.basename, "_gene_fpkm.tsv.gz", sep="")))
  write.tsv(isoform.fpkm,
    file=gzfile(paste(output.basename, "_isoform_fpkm.tsv.gz", sep="")))
}

system(paste("mkdir -p", output.path))

# r = get.cufflinks.results(
#   paste(input.path, "/Murray_050912/cehm6", sep=""))

for(a in list.files(input.path)) {
  cat(a, "\n")
  convert.dir(paste(input.path, a, sep=""),
    paste(output.path, a, sep=""))
}


