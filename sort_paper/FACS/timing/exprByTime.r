# Plots expression from various experiments, of genes
# which are enriched at particular times.

source("git/utils.r")
source("git/data/name_convert.r")

# genes expressed at particular times
gene.expr.time =
  read.tsv("git/sort_paper/FACS/timing/geneExprTime.tsv")

sd.cutoff = quantile(gene.expr.time$sd, 0.25)
td1 = gene.expr.time[ gene.expr.time$sd <= sd.cutoff , ]

# read data; we recompute enrichments relative to singlets
r = read.tsv("git/cluster/readsPerMillion.tsv")
pseudocount = 3
r = log2( as.matrix(r) + pseudocount )

# the Spencer data
load("data/expression/spencer.expr.Rdata")
s1 = spencer.expr[[1]]
colnames(s1)[[7]] = "LE.reference"
rownames(s1) = s1$name
s1 = as.matrix(s1[ , 4:33 ])

# averaged singlet control
singlet.average = (r[,"cnd-1 singlets"] + r[,"pha-4 singlets"]) / 2

# Plots profile of enrichment.
plot.time.profile = function(a, b, main, ylab="enrichment") {
  write.status(main)
  plot(td1[,"mean"],
    a[ rownames(td1) ] - b[ rownames(td1) ],
    main = main, xlab = "time", ylab=ylab,
    pch=20, cex=0.5, col="#00000080")
  abline(h = 0, col="#606060")
}

if (TRUE) {
  pdf("git/sort_paper/FACS/timing/exprByTime.pdf",
    width=7.5, height=10)
  par(mfrow=c(4,2))

  samples = grep("\\(\\+\\)", colnames(r), value=TRUE)
  samples = sort(grep("ceh-6.*hlh-16", samples, value=TRUE, invert=TRUE))

  for(s in samples) {

    plot.time.profile(r[,s], singlet.average,
      paste(s, "vs singlets"))
  }

  # Spencer data
  plot.time.profile(s1[,"EE.GLP"], s1[,"EE.reference"], "EE.GLP")
  plot.time.profile(s1[,"EE.BAG"], s1[,"EE.reference"], "EE.BAG")

  plot.time.profile(s1[,"LE.panneural"], s1[,"LE.reference"],
    "LE.panneural")
  plot.time.profile(s1[,"LE.AVA"], s1[,"LE.reference"],
    "LE AVA neurons")
  plot.time.profile(s1[,"LE.AVE"], s1[,"LE.reference"],
    "LE AVE neurons")
  plot.time.profile(s1[,"LE.A.class"], s1[,"LE.reference"],
    "LE A class neurons")
  plot.time.profile(s1[,"LE.bwm"], s1[,"LE.reference"],
    "LE body wall muscle")
  plot.time.profile(s1[,"LE.coelomocytes"], s1[,"LE.reference"],
    "LE coelomocytes")
  plot.time.profile(s1[,"LE.dopaminergic"], s1[,"LE.reference"],
    "LE dopaminergic neurons")
  plot.time.profile(s1[,"LE.hypodermis"], s1[,"LE.reference"],
    "LE hypodermis")
  plot.time.profile(s1[,"LE.AVE"], s1[,"LE.intestine"],
    "LE intestine")
  plot.time.profile(s1[,"LE.PhM"], s1[,"LE.reference"],
    "LE pharyngeal muscle")

  dev.off()
}




