# Summarizes which clusters have anatomy annotation (either
# anatomy terms, or genes in an anatomy-specific cluster.)

source("git/utils.r")

anatomy.term = read.tsv(
  "git/sort_paper/enrichment/summary/anatomyEnrichment/hier.300.clusters.tsv")
wb.cluster = read.tsv(
  "git/sort_paper/enrichment/summary/wormbaseCluster/hier.300.clusters.tsv")

# pick WormBase clusters which are anatomy-specific (e.g. from Spencer data)
wb.cluster = wb.cluster[ grep("WBPaper00037950", wb.cluster$group) , ]

# embryonic subset of clusters
wb.cluster.embryonic = wb.cluster[ (grepl("Embryo|embryo", wb.cluster$group) | grepl("Embryo|embryo", wb.cluster$group.name)) , ]

# remove some trivial annotation
anatomy.term = anatomy.term[ anatomy.term$group.name != "Tissue" , ]

# only keep one term per cluster
anatomy.term = anatomy.term[ ! duplicated(anatomy.term$cluster) , ]
wb.cluster = wb.cluster[ ! duplicated(wb.cluster$cluster) , ]

# preferentially use anatomy term annotation
wb.cluster = wb.cluster[ ! ( wb.cluster$cluster %in% anatomy.term$cluster ) , ]

# compact name
anatomy.term$name.short = anatomy.term$group.name
wb.cluster$name.short = wb.cluster$group

# print some stats
cat(paste("num. clusters with anatomy annotation =",
  nrow(anatomy.term), "\n"))
cat(paste("num. add'l with anatomy-specific expression clusters =",
  nrow(wb.cluster), "\n"))

anatomy.annotated = rbind(anatomy.term, wb.cluster)
anatomy.annotated = anatomy.annotated[ order(anatomy.annotated$p.corr) , ]

write.tsv(anatomy.annotated,
  "git/sort_paper/enrichment/clustersWithAnatomyAnnotation.tsv")

wb.cluster.nonembryonic =
  wb.cluster[ ! ( wb.cluster$cluster %in% wb.cluster.embryonic$cluster ) , ]

