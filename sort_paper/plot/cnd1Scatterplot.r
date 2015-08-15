# Plots cnd-1 enrichments.

source("git/plot/label_panel.r")
source("git/utils.r")
source("git/data/name_convert.r")

# expression data
rpm = read.tsv("git/cluster/readsPerMillion.tsv")

# read data, just including the FACS data, and only
# including cases in which all the data is present
r = read.tsv("git/cluster/readRatios.tsv")
r = r[,c(1:23)]
r = r[ apply(!is.na(r), 1, all) , ]

rpm.facs = rpm[ , c(1:23,38:54,65:68) ]

# omitting ribosomal RNA
rpm.facs = rpm.facs[ rownames(rpm.facs) != "ribosomal_RNA" , ]

# Plots a scatterplot of two different samples.
plot.scatter = function(s1, s2) {

  lim = range(c(r[, c(s1, s2)]))
  plot(r[ , s1], r[ , s2],
    col="#00000040", cex=0.5,
    pch=183, font=5, xlim=lim, ylim=lim, xaxt="n", yaxt="n",
    xlab=expression(italic("cnd-1") * " enrichment, replicate 1"),
      ylab=expression(italic("cnd-1") * " enrichment, replicate 2"))
  axis(1)
  axis(2)
  abline(0,1, col="#00000060")
}

pdf("git/sort_paper/plot/cnd1Scatterplot.pdf",
  width=4, height=4)   # was 5x5
# png("git/sort_paper/plot/cnd1Scatterplot.png",
#   width=4, height=4, units="in", res=300)

par(mar=c(4,4,0.5,0.5)+0.1)

plot.scatter("cnd-1 8/19", "cnd-1 12/14")
# plot.scatter("cnd-1 8/19", "cnd-1 1/4")
# plot.scatter("cnd-1 12/14", "cnd-1 1/4")

dev.off()

