# Plots expression from various experiments, of genes
# which are enriched at particular times.

source("git/utils.r")
source("git/data/name_convert.r")

source("git/sort_paper/FACS/timing/proportionOnFromImaging.r")

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
plot.time.profile = function(a, b, main, ylab="Enrichment") {
  write.status(main)
  plot(td1[,"mean"],
    a[ rownames(td1) ] - b[ rownames(td1) ],
    main = main, xlab = "Time (minutes)", ylab=ylab,
    ylim=c(-5,5), xaxt="n", yaxt="n",
    pch=183, font=5, cex=1, col="#00000080")
  axis(1)
  axis(2)
  abline(h = 0, col="#606060")
}

# Plots expression versus time for all the sorting experiments,
# and experiments from Spencer et al.
plot.all.expr.by.time = function() {
  pdf("git/sort_paper/FACS/timing/exprByTime.pdf",
    width=7.5, height=10)
  par(mfrow=c(4,2))

  samples = grep("\\(\\+\\)", colnames(r), value=TRUE)
  samples = sort(grep("ceh-6.*hlh-16", samples, value=TRUE, invert=TRUE))

  for(s in samples) {
    plot.time.profile(r[,s], singlet.average,
      paste(s, "vs singlets"))
  }

  # positive vs. negative
  for(s in samples) {
    neg = sub("\\(\\+\\)", "\\(\\-\\)", s)
    if (neg %in% colnames(r)) {
      plot.time.profile(r[,s], r[,neg], paste(s, "vs", neg))
    }
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

# simple smoothing
sm = function(p) {
  m = loess(proportion.on ~ t, p, alpha = 0.1)
  p$proportion.on = predict(m, p)
  p$proportion.on[ p$proportion.on < 0 ] = 0
  p
}

# Plots the expression profile, and how many cells were expressing
# at particular times, for a few genes.
#plot.expr.and.percent.on = function() {

on.off.c = "#0000ffb0"

  # cnd-1
  pdf("git/sort_paper/FACS/timing/exprByTime_cnd1.pdf",
    width=8, height=4.5)
  par(mar=c(5,4,4,4) + 0.1)

  plot.time.profile(r[,"cnd-1 8/19 (+)"], singlet.average, "")

  p = sm(get.proportion.on("SCD20110731_RW10434_L1_smooth", 50, 1000))
  par(new=TRUE)
  ymax=signif(max(p$proportion.on), 1)
  p = p[p$t <= 400,]   # only show part including data
  plot(p$t, p$proportion.on, type="l", lwd=3, col=on.off.c,
    xlim=range(td1$mean), ylim=c(-ymax, ymax), xaxt="n", yaxt="n",
    xlab="", ylab="")
  axis(4, at = c(0, 0.5, 1) * ymax, col=on.off.c, col.axis=on.off.c)
  mtext("Proportion of expressing cells", col=on.off.c, side=4, line=2)
  dev.off()

  # pros-1
  pdf("git/sort_paper/FACS/timing/exprByTime_pros1.pdf",
    width=8, height=4.5)
  par(mar=c(5,4,4,4) + 0.1)

  plot.time.profile(r[,"pros-1 (+)"], singlet.average, "")

  p = sm(get.proportion.on("SCD20120426_JIM122_L2_smooth", 50, 1000))
  par(new=TRUE)
  ymax=signif(max(p$proportion.on), 1)
  p = p[p$t <= 400,]   # XXX only show part including data
  plot(p$t, p$proportion.on, type="l", lwd=3, col=on.off.c,
    xlim=range(td1$mean), ylim=c(-ymax, ymax), xaxt="n", yaxt="n",
    xlab="", ylab="")
  axis(4, at = c(0, 0.5, 1) * ymax, col=on.off.c, col.axis=on.off.c)
  mtext("Proportion of expressing cells", col=on.off.c, side=4, line=2)
  dev.off()

if (FALSE) {
  pdf("git/sort_paper/FACS/timing/proportionOn.pdf",
    width=8, height=4.5)
  par(mar=c(5,4,4,4) + 0.1)

  p = sm(get.proportion.on("SCD20110726_RW10713_L1_smooth", 50, 1000))
  plot(p$t, p$proportion.on, type="l")

  dev.off()
}

#}

# plot.all.expr.by.time()

