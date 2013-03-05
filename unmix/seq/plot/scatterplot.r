# Does some scatterplots of sorted fractions.

output.path = "git/unmix/seq/plot/scatterplot/"

system(paste("mkdir -p ", output.path))

options(stringsAsFactors = FALSE)

experimentNames = read.table("git/unmix/seq/quant/experimentNames.tsv",
  sep="\t", header=TRUE, row.names=1)

count.path = "git/unmix/seq/quant/readsPerMillion"

# get the reads from all experiments
r1 = as.matrix(
  read.table(paste(count.path, "readsPerMillion_Murray_050912.tsv", sep="/"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))
r2 = as.matrix(
  read.table(paste(count.path, "readsPerMillion_092812.tsv", sep="/"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))
r3 = as.matrix(
  read.table(paste(count.path, "readsPerMillion_20110922.tsv", sep="/"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

r = log2(1 + cbind(r1, r2, r3))    # ??? add something smaller than 1 here?
r = r[ , order(colnames(r)) ]
colnames(r) = experimentNames[colnames(r),"name"]

r.ds = r[,c("ceh-6 (+) hlh-16 (+)", "ceh-6 (+) hlh-16 (-)",
  "ceh-6 (-) hlh-16 (+)", "ceh-6 (-) hlh-16 (-)")]
r.ds = cbind(r.ds, "mean of double-sorted samples" = apply(r.ds, 1, mean))


# Draws a scatterplot.
draw.scatter.1 = function(r, x.name, y.name, include.correlation=FALSE) {
  x = r[,x.name]
  y = r[,y.name]
  lim = c(0, max(x, y))
#  main = if (include.correlation)
#    paste("cor =", round(cor(x,y), 3))
#  else ""
  main = expression("expression, log"[2](1+coverage))

  plot(r[,x.name], r[,y.name], xlab=x.name, ylab=y.name,
    main = main, xlim=lim, ylim=lim, pch=20, col="#00000080",
    cex=1.5, cex.axis=1.5, cex.lab=2, cex.main=1.6)
  abline(0,1, lwd=3, col="#00000060")
}

# Draws a scatterplot; variation on above.
draw.scatter.2 = function(x, y, x.name, y.name, include.correlation=FALSE) {
  lim = c(min(x,y), max(x, y))
  main = expression("log"[2](enrichment))

  plot(x, y, xlab=x.name, ylab=y.name,
    main = main, xlim=lim, ylim=lim, pch=20, col="#00000080",
    cex=1.5, cex.axis=1.5, cex.lab=2, cex.main=1.6)
  abline(0,1, lwd=3, col="#00000060")
}

# Draws a scatterplot, coloring genes which are far away from
# a threshold.
draw.scatter.diag.threshold =
    function(x, y, x.name, y.name, cutoff = 2) {
  lim = c(min(x,y), max(x, y))
  main = expression("log"[2](enrichment))

  colors = ifelse(y-x > cutoff, "#80000080",
    ifelse(y-x < -cutoff, "#00008080", "#00000080"))

  plot(x, y, xlab=x.name, ylab=y.name,
    main = main, xlim=lim, ylim=lim, pch=20, col=colors,
    cex=1.5, cex.axis=1.5, cex.lab=2, cex.main=1.6)
  abline(0,1, lwd=3, col="#00000060")
  abline(2,1, lwd=3, col="#80000060")
  abline(-2,1, lwd=3, col="#00008060")
}

scatterplots1 = function() {
  draw.scatter = function(x,y) draw.scatter.1(r,x,y)

  png(file=paste(output.path, "/", "plusMinus.png", sep=""),
    width=900, height=600)
  par(mar=c(5,5,4,2)+0.1)
  par(mfrow=c(2,3))
  draw.scatter("ceh-6 (+)", "ceh-6 (-)")
  draw.scatter("ceh-26 (+)", "ceh-26 (-)")
  draw.scatter("ceh-27 (+)", "ceh-27 (-)")

  draw.scatter("ceh-36 (+)", "ceh-36 (-)")
  draw.scatter("F21D5.9 (+)", "F21D5.9 (-)")
  draw.scatter("mls-2 (+)", "mls-2 (-)")
  dev.off()

  png(file=paste(output.path, "/", "cnd-1_reps.png", sep=""),
    width=900, height=600)
  par(mar=c(5,5,4,2)+0.1)
  par(mfrow=c(2,3))
  draw.scatter("cnd-1 8/19 (+)", "cnd-1 12/14 (+)")
  draw.scatter("cnd-1 8/19 (+)", "cnd-1 1/4 (+)")
  draw.scatter("cnd-1 12/14 (+)", "cnd-1 1/4 (+)")

  draw.scatter("cnd-1 8/19 (-)", "cnd-1 12/14 (-)")
  draw.scatter("cnd-1 8/19 (-)", "cnd-1 1/4 (-)")
  draw.scatter("cnd-1 12/14 (-)", "cnd-1 1/4 (-)")
  dev.off()

  png(file=paste(output.path, "/", "pha-4_reps.png", sep=""),
    width=900, height=600)
  par(mar=c(5,5,4,2)+0.1)
  par(mfrow=c(2,3))
  draw.scatter("pha-4 5/9 (+)", "pha-4 9/1 (+)")
  draw.scatter("pha-4 5/9 (+)", "pha-4 12/9 (+)")
  draw.scatter("pha-4 9/1 (+)", "pha-4 12/9 (+)")

  draw.scatter("pha-4 5/9 (-)", "pha-4 9/1 (-)")
  draw.scatter("pha-4 5/9 (-)", "pha-4 12/9 (-)")
  draw.scatter("pha-4 9/1 (-)", "pha-4 12/9 (-)")
  dev.off()

  png(file=paste(output.path, "/", "replicateScatter.png", sep=""),
    width=900, height=600)
  par(mar=c(5,5,4,2)+0.1)
  par(mfrow=c(2,3))
  draw.scatter("cnd-1 8/19 (-)", "cnd-1 8/19 (+)")
  draw.scatter("cnd-1 12/14 (-)", "cnd-1 12/14 (+)")
  draw.scatter("cnd-1 1/4 (-)", "cnd-1 1/4 (+)")

  draw.scatter("pha-4 5/9 (-)", "pha-4 5/9 (+)")
  draw.scatter("pha-4 9/1 (-)", "pha-4 9/1 (+)")
  draw.scatter("pha-4 12/9 (-)", "pha-4 12/9 (+)")
  dev.off()

  png(file=paste(output.path, "/", "ceh-6_hlh-16.png", sep=""),
    width=900, height=600)
  par(mar=c(5,5,4,2)+0.1)
  par(mfrow=c(2,3))
  draw.scatter("ceh-6 (+)", "ceh-6 (+) hlh-16 (+)")
  draw.scatter("ceh-6 (+)", "ceh-6 (+) hlh-16 (-)")
  draw.scatter("ceh-6 (+)", "ceh-6 (-) hlh-16 (+)")

  draw.scatter("ceh-6 (+)", "ceh-6 (-) hlh-16 (-)")
  draw.scatter("ceh-6 (-)", "ceh-6 (-) hlh-16 (+)")
  draw.scatter("ceh-6 (-)", "ceh-6 (-) hlh-16 (-)")

  dev.off()

  png(file=paste(output.path, "/", "ceh-6_hlh-16_2.png", sep=""),
    width=800, height=800)
  par(mar=c(5,5,4,2)+0.1)
  par(mfrow=c(2,2))

  draw.scatter("hlh-16 (+)", "ceh-6 (-) hlh-16 (+)")
  draw.scatter("hlh-16 (+)", "ceh-6 (-) hlh-16 (+)")

  draw.scatter("ceh-6 (+) hlh-16 (+)", "ceh-6 (-) hlh-16 (-)")
  draw.scatter("ceh-6 (+) hlh-16 (-)", "ceh-6 (-) hlh-16 (+)")
  dev.off()
}

scatterplots2 = function() {

  for(g in c("ceh-26", "ceh-6", "pha-4 9/1", "ttx-3")) {
    g1 = sub("/", "_", g)
    png(file=paste(output.path, "/", g1, ".png", sep=""),
      width=800, height=800)
    par(mar=c(5,5,4,2)+0.1)
    draw.scatter.1(r, paste(g, "(-)"), paste(g, "(+)"))
    dev.off()
  }

  png(file=paste(output.path, "/", "ceh-6_hlh-16.png", sep=""),
    width=800, height=800)
  par(mar=c(5,5,4,2)+0.1)
  draw.scatter.1(r.ds, "mean of double-sorted samples",
    "ceh-6 (+) hlh-16 (+)")
  dev.off()
}

# Plots relating to replicates.
replicate.scatterplots = function() {

  png(file=paste(output.path, "/pha4_replicates.png", sep=""),
    width=1200, height=400)
  par(mfrow=c(1,3))
  par(mar=c(5,5,4,2)+0.1)

  draw.scatter.2(r[,"pha-4 5/9 (+)"] - r[,"pha-4 5/9 (-)"],
    r[,"pha-4 9/1 (+)"] - r[,"pha-4 9/1 (-)"],
    "pha-4 (5/9)", "pha-4 (9/1)", "")
  draw.scatter.2(r[,"pha-4 5/9 (+)"] - r[,"pha-4 5/9 (-)"],
    r[,"pha-4 12/9 (+)"] - r[,"pha-4 12/9 (-)"],
    "pha-4 (5/9)", "pha-4 (12/9)", "")
  draw.scatter.2(r[,"pha-4 9/1 (+)"] - r[,"pha-4 9/1 (-)"],
    r[,"pha-4 12/9 (+)"] - r[,"pha-4 12/9 (-)"],
    "pha-4 (9/1)", "pha-4 (12/9)", "")

  dev.off()

  png(file=paste(output.path, "/cnd1_replicates.png", sep=""),
    width=1200, height=400)
  par(mfrow=c(1,3))
  par(mar=c(5,5,4,2)+0.1)

  draw.scatter.2(r[,"cnd-1 8/19 (+)"] - r[,"cnd-1 8/19 (-)"],
    r[,"cnd-1 12/14 (+)"] - r[,"cnd-1 12/14 (-)"],
    "cnd-1 (8/19)", "cnd-1 (12/14)", "")
  draw.scatter.2(r[,"cnd-1 8/19 (+)"] - r[,"cnd-1 8/19 (-)"],
    r[,"cnd-1 1/4 (+)"] - r[,"cnd-1 1/4 (-)"],
    "cnd-1 (8/19)", "cnd-1 (1/4)", "")
  draw.scatter.2(r[,"cnd-1 12/14 (+)"] - r[,"cnd-1 12/14 (-)"],
    r[,"cnd-1 1/4 (+)"] - r[,"cnd-1 1/4 (-)"],
    "cnd-1 (12/14)", "cnd-1 (1/4)", "")

  dev.off()
}

# Scatterplots of genes enriched / depleted in particular fractions.
enriched.depleted.scatterplots = function() {

  png(file=paste(output.path, "/pha4_enrichment.png", sep=""),
    width=800, height=800)
  par(mar=c(5,5,4,2)+0.1)
  draw.scatter.diag.threshold(r[,"pha-4 9/1 (-)"], r[,"pha-4 9/1 (+)"],
    "pha-4 9/1 (-)", "pha-4 9/1 (+)", 2)
  dev.off()
}

# scatterplots1()
# scatterplots2()
# replicate.scatterplots()

enriched.depleted.scatterplots()

