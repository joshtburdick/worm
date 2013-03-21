# Compares correlation with a TF, and motifs of that TF.

source("git/utils.r")

x = read.tsv("git/unmix/seq/cluster/readsFACSandTS.tsv")

# remove cases with no variance
x = x[ apply(x,1,var) > 0 , ]

# motif counts
motif = t(as.matrix(read.tsv(
  gzfile("git/tf/motif/motifCount/motifs_5kbUpstream.tsv.gz"))))

# only keep genes which are in both lists
g = intersect(rownames(x), rownames(motif))
x = x[g,]
motif = motif[g,]

# instances of orthologs
ortho = read.table("git/tf/wb.ortholog.tsv", header=TRUE, as.is=TRUE)

ortho = ortho[ ortho$public_name %in% rownames(x), ]

# just keep cases in which there's a motif, and only keep
# the motif from the "best" homolog (by blastp evalue)
ortho = ortho[ ortho$motif != "" , ]
ortho = ortho[
  order(ortho$sequence_name,
    ortho$homolog_species,
    ortho$homolog_blastp_evalue) , ]
ortho = ortho[ !duplicated(ortho[ , c("sequence_name", "homolog_species"),]) , ]

# XXX some of these motifs aren't included (I'm not sure why),
# so omit those cases
ortho = ortho[ ortho$motif %in% colnames(motif) , ]


tfs = intersect(rownames(x), ortho$public_name)

tf.cor = cor(t(x), t(x[tfs,]))

# also omit cases in which expression is missing (presumably
# because of a naming issue)
ortho = ortho[ ortho$public_name %in% colnames(tf.cor) , ]

cor.and.motif = cor( tf.cor[,ortho$public_name], motif[,ortho$motif] )


plot.it.1 = function(i) {
  plot(motif[,ortho[i,"motif"] ], tf.cor[,ortho[i,"public_name"] ])
}




plot.it.2 = function(gene.name, motif.name, n) {

  # which genes have high and low correlation, respectively
  hi = order(tf.cor[,gene.name ], decreasing=TRUE)[1:n]
  lo = order(abs(tf.cor[,gene.name ]), decreasing=FALSE)[1:n]
#  lo = setdiff(c(1:nrow(tf.cor)), hi)
  hi.motif.count = sqrt(motif[hi,motif.name ])
  lo.motif.count = sqrt(motif[lo,motif.name ])
  xlim = c(0, max(lo.motif.count, hi.motif.count))

  par(mfrow=c(2,1))
  hist(hi.motif.count, xlim=xlim, col="grey")
  hist(lo.motif.count, xlim=xlim, col="grey")
  s = (t.test(hi.motif.count, lo.motif.count))
  print(s)
}


plot.it.2("daf-19", "RFX2_DBD", 1000)
plot.it.2("ceh-36", "OTX1_DBD", 1000)
plot.it.2("pha-4", "MA0148.1", 1000)

