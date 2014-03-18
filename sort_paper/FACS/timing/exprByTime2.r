# Plots expression from various experiments, of genes
# which are enriched at particular times.
# This version groups genes at particular times, rather
# than plotting a timecourse.

library("vioplot")

source("git/utils.r")
source("git/data/name_convert.r")

# read data; we recompute enrichments relative to singlets
r = read.tsv("git/cluster/readsPerMillion.tsv")
pseudocount = 3
r = log2( as.matrix(r) + pseudocount )

# genes expressed at particular times
gene.expr.time =
  read.tsv("git/sort_paper/FACS/timing/geneExprTime.tsv")

sd.cutoff = quantile(gene.expr.time$sd, 0.25)
td1 = gene.expr.time[ gene.expr.time$sd <= sd.cutoff , ]
td1 = td1[ rownames(td1) %in% rownames(r) , ]

# define "categories" of times
time.name = list(
  Early = rownames(td1)[ td1$mean < 300 ],
  Middle = rownames(td1)[ td1$mean >= 300 & td1$mean < 450 ],
  Late = rownames(td1)[ td1$mean >= 450 ])


# averaged singlet control
singlet.average = (r[,"cnd-1 singlets"] + r[,"pha-4 singlets"]) / 2

# Plots profile of enrichment, as violin plots.
plot.time.profile = function(a, b, main, ylab="enrichment") {
  write.status(main)

  x = a[ rownames(td1) ] - b[ rownames(td1) ]
# browser()
  par(mar=c(5,4,4,2)-0.5)
  vioplot(x[ time.name[["Early"]] ],
    x[ time.name[["Middle"]] ],
    x[ time.name[["Late"]] ],
    names = c("Early", "Middle", "Late"),
    col="grey", ylim=c(-6,6))
  title(main)
}

if (TRUE) {
  pdf("git/sort_paper/FACS/timing/exprByTime2.pdf",
    width=7.5, height=10)
  par(mfrow=c(4,3))

  samples = grep("\\(\\+\\)", colnames(r), value=TRUE)
  samples = sort(grep("ceh-6.*hlh-16", samples, value=TRUE, invert=TRUE))

  # positive vs. negative
  for(s in samples) {
    neg = sub("\\(\\+\\)", "\\(\\-\\)", s)
    if (neg %in% colnames(r)) {
      plot.time.profile(r[,s], r[,neg], 
      paste(s, "enrichment"))
    }
  }

  dev.off()
}

