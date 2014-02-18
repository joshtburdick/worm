# Measures enrichment of sort markers in
# FACS-sorted samples.

source("git/plot/label_panel.r")
source("git/utils.r")

# expression data
rpm = read.tsv("git/cluster/readsPerMillion.tsv")

r = read.tsv("git/cluster/readRatios.tsv")

r.sort.only = r[,c(1:23)]
rpm.facs = rpm[ , c(1:23,38:54,65:68) ]

# omitting ribosomal RNA
rpm.facs = rpm.facs[ rownames(rpm.facs) != "ribosomal_RNA" , ]

# used for matching up gene names
gene.ids = read.csv("data/wormbase/geneIDs.WS224.csv.gz",
  header=FALSE, as.is=TRUE)

# Matches up gene names.
match.names = function(gene.names) {
  g = gene.ids[ gene.ids[,3] %in% gene.names , ]
  intersect(rownames(r), union(g[,2], g[,3]))
}

get.sort.marker.enrichment = function() {
  sample.names = colnames(r)[1:18]
  # ??? possibly omit hlh-16?
  sample.names = sample.names[ !(sample.names %in% c("mir-57")) ]
  genes = gsub(" [^ ]+", "", sample.names)
  enrich = r[ cbind(genes, sample.names) ]  
  names(enrich) = sample.names
  enrich[ !is.na(enrich) ]
}

# which genes are promoter fusions
promoter.fusion = c("ceh-27", "ceh-36", "F21D5.9", "hlh-16", "irx-1")

pdf("git/sort_paper/plot/enrichmentInSortFractions.pdf",
  width=7, height=5)

sort.marker.enrichment = sort(get.sort.marker.enrichment())


par(mar=c(7,4,4,2)+0.1)
barplot(sort.marker.enrichment, space=0.5, las=2,
  main="Enrichment of sort marker in sorted fraction",
  ylab="Enrichment",
  col=ifelse(names(sort.marker.enrichment) %in% promoter.fusion,
    "#c00000", "#80ff80"))
legend("topleft", c("promoter fusion", "protein fusion"),
  fill=c("#c00000", "#80ff80"))
abline(h=0)

dev.off()
