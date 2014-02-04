# Counts of the number of genes enriched in
# particular FACS-sorted samples.

source("git/plot/label_panel.r")
source("git/utils.r")
source("git/data/name_convert.r")

# expression data
rpm = read.tsv("git/cluster/readsPerMillion.tsv")

r = read.tsv("git/cluster/readRatios.tsv")

r.sort.only = r[,c(1:23)]
rpm.facs = rpm[ , c(1:23,38:54,65:68) ]

# omitting ribosomal RNA
rpm.facs = rpm.facs[ rownames(rpm.facs) != "ribosomal_RNA" , ]

# pharyngeal genes (Gaudet 2003)
pha4.early = read.csv("data/expression/gaudet/Ph-E_genes.csv",
  as.is=TRUE)
pha4.late = read.csv("data/expression/gaudet/Ph-L_genes.csv",
  as.is=TRUE)
rownames(pha4.early) = pha4.early$gene
rownames(pha4.late) = pha4.late$gene

pha4.early = rename.gene.names(pha4.early)
pha4.late = rename.gene.names(pha4.late)

pha4.replicates = c("pha-4 9/1", "pha-4 5/9", "pha-4 12/9")

# Plots a histogram of genes.
plot.it = function(r) {
  r = as.matrix(r)
  r = r[ order(apply(r,1,mean)), ]
  barplot(t(r), beside=TRUE, las=3, cex.names=0.8, yaxt="n",
    ylim=c(-2,6))
  axis(side=2)
  mtext("Enrichment", cex=1.1, side=2, line=3)
  abline(h=0, col="darkgrey")
}

pdf("git/sort_paper/plot/pha4Enrichment.pdf",
  width=10, height=7)
par(mfcol=c(2,1))

if (FALSE) {
hist(r[ rownames(pha4.early) , "pha-4 9/1" ],
  main="Early pharyngeal genes",
  xlab="enrichment", xlim=c(-2,6), col="grey")
label.panel("a)")
hist(r[ rownames(pha4.late) , "pha-4 9/1" ],
  main="Late pharyngeal genes",
  xlab="enrichment", xlim=c(-2,6), col="grey")
label.panel("b)")
}

plot.it(r[ rownames(pha4.early), pha4.replicates ])
mtext("Early pharyngeal genes", cex=1.3)
label.panel("a)")
plot.it(r[ rownames(pha4.late), pha4.replicates ])
mtext("Late pharyngeal genes", cex=1.3)
label.panel("b)")
dev.off()

