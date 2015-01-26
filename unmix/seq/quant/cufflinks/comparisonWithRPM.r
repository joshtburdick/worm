# Compares FPKM and RPM.

source("git/utils.r")

output.dir = "git/unmix/seq/quant/cufflinks/comparisonWithRPM/"

# read ratios, including intermediate results
source("git/cluster/readRatios.r")
rpm.enrich = r.ft

cuff.result = "git/unmix/seq/quant/cufflinks/cufflinks_summary_20140121/"

gene.fpkm = cbind(
  read.tsv(paste0(cuff.result, "Murray_050912_gene_fpkm.tsv.gz")),
  read.tsv(paste0(cuff.result, "Murray_092812_gene_fpkm.tsv.gz")),
  read.tsv(paste0(cuff.result, "Murray_20110922_gene_fpkm.tsv.gz")))

g = intersect(rownames(readsPerMillion), rownames(gene.fpkm))
rpm.enrich = rpm.enrich[g,]
gene.fpkm = gene.fpkm[g,]
readsPerMillion = readsPerMillion[g,]

# one way of computing enrichment
enrich = function(pos, neg, eps=1) {
  r = log2(eps + pos) - log2(eps + neg)
  r
}

# Plots a comparison of enrichment, and writes out a table
# of very-different enrichment ratios (and the corresponding
# expression numbers which went into them.)
enrich.comparison = function(gene) {
  pos.name = paste(gene, "(+)")
  neg.name = paste(gene, "(-)")
  main = gsub("/", "_", gene)

  pos.fpkm.name = rownames(experimentNames)[experimentNames$name==pos.name]
  neg.fpkm.name = rownames(experimentNames)[experimentNames$name==neg.name]
# cat("names are ", pos.fpkm.name, neg.fpkm.name, "\n")

  r = data.frame(cbind(rpm.pos = readsPerMillion[ , pos.name ],
    rpm.neg = readsPerMillion[ , neg.name ],
    fpkm.pos = gene.fpkm[ , pos.fpkm.name ],
    fpkm.neg = gene.fpkm[ , neg.fpkm.name ]))

  r$rpm.enrich = enrich(r$rpm.pos, r$rpm.neg)
  r$fpkm.enrich = enrich(r$fpkm.pos, r$fpkm.neg)
  r$enrich.diff = r$fpkm.enrich - r$rpm.enrich
  r$enrich.diff[ is.na(r$enrich.diff) ] = 0

  png(paste0(output.dir, "/", main, ".png"), width=1000, height=1000)
  plot(r$rpm.enrich, r$fpkm.enrich,
    type="p", pch=20, col="#00000060", main=gene)
  dev.off()

  r1 = r[ abs(r$enrich.diff) >= 1 , ]
  r1 = r1[ order(r1$enrich.diff), ]
  write.tsv(r1, paste0(output.dir, "/", main, ".tsv"))
}

system(paste0("mkdir -p ", output.dir))

for(a in c("ceh-26", "ceh-27", "ceh-36", "ceh-6", "mls-2", "pal-1",
    "pha-4 5/9", "pha-4 9/1", "pha-4 12/9",
    "cnd-1 1/4", "cnd-1 8/19", "cnd-1 12/14",
    "F21D5.9", "mir-57", "unc-130", "ttx-3"))
  enrich.comparison(a)



