# Summarizes which clusters have anatomy annotation (either
# anatomy terms, or genes in an anatomy-specific cluster.)

source("git/utils.r")

anatomy.term = read.tsv(
  "git/sort_paper/enrichment/summary/anatomyEnrichment/hier.300.clusters.tsv")
wb.cluster.all = read.tsv(
  "git/sort_paper/enrichment/summary/wormbaseCluster/hier.300.clusters.tsv")

# pick WormBase clusters which are anatomy-specific (e.g. from Spencer data)
# XXX rename this?
wb.cluster = wb.cluster.all[ grep("WBPaper00037950", wb.cluster.all$group) , ]

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

##### Finding clusters which aren't tissue specific
# first, omit anything annotated already
wb.cluster.non.ts =
  wb.cluster.all[ !(wb.cluster.all$cluster %in% anatomy.term$cluster) , ]

# remove tissue-specific experiments
ts.pattern = "(WBPaper(00037950|00031003|00030839|00026980|00025141|00036375|00031532|00024671|00025032|00029359))|cgc5428|cgc6390"
wb.cluster.non.ts =
  wb.cluster.non.ts[ ! grepl(ts.pattern, wb.cluster.non.ts$group) , ]
# the procedure is to hack the above pattern, until there's nothing
# tissue-specific left in this file
write.tsv(wb.cluster.non.ts, "git/sort_paper/enrichment/WBclusterNonTS.tsv")

cat("num. non-ts clusters =", length(unique(wb.cluster.non.ts$cluster)), "\n")

