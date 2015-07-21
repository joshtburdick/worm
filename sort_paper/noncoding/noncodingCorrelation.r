# Correlation of non-coding genes with nearby genes.

source("git/utils.r")
source("git/plot/label_panel.r")

# the read data, as enrichment ratios
r = read.tsv("git/cluster/readRatios.tsv")
r.sort.only = r[,c(1:23)]

# locations of genes
gene.loc = read.table("git/data/seq/merged_genes_WS220.bed",
  sep="\t", as.is=TRUE)
colnames(gene.loc) = c("chr", "a", "b", "gene", "score", "strand")

# names of noncoding genes
noncoding = grep("^(anr|linc)-", rownames(r.sort.only), value=TRUE)

# alternative way of matching gene names, using WormBase annotation
gene.ids = read.csv(gzfile("data/wormbase/geneIDs.WS224.csv.gz"),
  as.is=TRUE, header=FALSE)
colnames(gene.ids) = c("wb.gene", "gene", "transcript")
gene.ids = gene.ids[ gene.ids[,"gene"] != "" , ]

# Utility to rename a vector of strings (e.g. gene names.)
# Args:
#   from, to - everything in "from" will be renamed "to"
#   a - the strings in question
# Returns: a, renamed
rename.gene.names = function(from, to, r) {
  i = intersect(r, from)
  r[ match(i, r) ] =
    to[ match(i, from) ]

  r
}

if (FALSE) {
# cases in which one gene is near a noncoding gene
nc.near.gene.1 = read.table("git/sort_paper/noncoding/noncodingNearGene.tsv",
  sep="\t", as.is=TRUE)

nc.near.gene = data.frame(nc.gene = nc.near.gene.1[,4], gene = nc.near.gene.1[,10],
  stringsAsFactors=FALSE)
nc.near.gene = nc.near.gene[ nc.near.gene$gene %in% rownames(r.sort.only) , ]
nc.near.gene = nc.near.gene[ nc.near.gene[,1] != nc.near.gene[,2] , ]

# nc.near.gene$gene = rename.gene.names(gene.ids[,"transcript"], gene.ids[,"gene"],
#   nc.near.gene$gene)

# find correlations

nc.near.gene$cor = NA
for(i in 1:nrow(nc.near.gene))
  nc.near.gene[i,"cor"] =
    cor(t(r.sort.only[nc.near.gene[i,"nc.gene"],]), t(r.sort.only[nc.near.gene[i,"gene"],]))
}

# correlations with closest gene
nc.closest.gene.1 = read.table("git/sort_paper/noncoding/closestToNoncoding.tsv",
  sep="\t", as.is=TRUE)
nc.closest.gene = data.frame(
  nc.gene = nc.closest.gene.1[,4],
  gene = nc.closest.gene.1[,10],
  dist = nc.closest.gene.1[,13], stringsAsFactors=FALSE)
nc.closest.gene =
  nc.closest.gene[ nc.closest.gene[,"gene"] %in% rownames(r.sort.only) , ]
nc.closest.gene = nc.closest.gene[
  nc.closest.gene[,"nc.gene"] != nc.closest.gene[,"gene"] , ]


nc.closest.gene = nc.closest.gene[ nc.closest.gene$dist <= 1000 , ]
nc.closest.gene$cor = NA
for(i in 1:nrow(nc.closest.gene)) {
  nc.closest.gene[i,"cor"] =
    cor(t(r.sort.only[nc.closest.gene[i,"nc.gene"],]), t(r.sort.only[nc.closest.gene[i,"gene"],]))
}

# correlations of all nearby genes
nearby.genes = read.table("git/sort_paper/noncoding/nearbyGenes.tsv",
  sep="\t", as.is=TRUE)
nearby.genes = nearby.genes[ , c(4,10) ]
colnames(nearby.genes) = c("a", "b")
nearby.genes = nearby.genes[ nearby.genes$a < nearby.genes$b , ]
nearby.genes$cor = NA
for(i in 1:nrow(nearby.genes)) {
  write.status(i)
  nearby.genes[i,"cor"] = cor(
    t(r.sort.only[nearby.genes[i,"a"],]),
    t(r.sort.only[nearby.genes[i,"b"],]))
}

anr.cor = nc.closest.gene[grep("anr-", nc.closest.gene$nc.gene),]
linc.cor = nc.closest.gene[grep("linc-", nc.closest.gene$nc.gene),]

pdf("git/sort_paper/noncoding/noncodingCorrelation.pdf", width=7, height=7)
par(mfrow=c(2,2))

hist(anr.cor$cor, breaks=20, xlim=c(-1,1), col="grey",
  main="anrRNA",
  xlab="Correlation with nearest gene")
mtext(paste("mean =", signif(mean(anr.cor$cor, na.rm=TRUE), 2)),
  cex=0.8)
label.panel("b)")

hist(linc.cor$cor, breaks=20, xlim=c(-1,1), col="grey",
  main="lincRNA",
  xlab="Correlation with nearest gene")
mtext(paste("mean =", signif(mean(linc.cor$cor, na.rm=TRUE), 2)),
  cex=0.8)
label.panel("c)")

hist(nearby.genes$cor, breaks=40, xlim=c(-1,1), col="grey",
  main="All adjacent genes",
  xlab="Correlation with nearest gene")
mtext(paste("mean =", signif(mean(nearby.genes$cor, na.rm=TRUE), 2)),
  cex=0.8)
label.panel("d)")

dev.off()

# cat("mean anr- correlation with nearest =", mean(anr.cor$cor, na.rm=TRUE),
#  "  (n=", sum(!is.na(anr.cor$cor)), ")\n")
# cat("mean linc- correlation with nearest =", mean(linc.cor$cor, na.rm=TRUE),
#   "  (n=", sum(!is.na(linc.cor$cor)), ")\n")


