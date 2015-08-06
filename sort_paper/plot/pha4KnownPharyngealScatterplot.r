# Plots enrichments of known pharyngeal genes,
# comparing replicates.

source("git/plot/label_panel.r")
source("git/utils.r")
source("git/data/name_convert.r")

# expression data
rpm = read.tsv("git/cluster/readsPerMillion.tsv")

r = read.tsv("git/cluster/readRatios.tsv")

r.sort.only = r[,c(1:23)]
rpm.facs = rpm[ , c(1:23,38:54,65:68) ]

# omitting ribosomal RNA
rpm.facs = rpm.facs[ rownames(rpm.facs) != "ribosomal_RNA" , ]

# pharyngeal genes (Gaudet 2003)
pha4.early = read.csv("data/expression/gaudet/Ph-E_genes.csv",
  as.is=TRUE)
pha4.late = read.csv("data/expression/gaudet/Ph-L_genes.csv",
  as.is=TRUE)
rownames(pha4.early) = pha4.early$gene
rownames(pha4.late) = pha4.late$gene

pha4.early = rename.gene.names(pha4.early)
pha4.late = rename.gene.names(pha4.late)

pha4.replicates = c("pha-4 5/9", "pha-4 9/1", "pha-4 12/9")
replicate.colors = hsv((0:2)/3, (3:1)/3, 1)

# Plots a scatterplot of two different samples.
plot.scatter = function(s1, s2) {

  g = c(rownames(pha4.early), rownames(pha4.late))
  lim = range(c(r[g, c(s1, s2)]))

  plot(r[ , s1], r[ , s2],
#    main = "Enrichment of known pharyngeal genes",
    col="#00000040", cex=0.5,
    pch=183, font=5, xlim=lim, ylim=lim, xaxt="n", yaxt="n",
    xlab=expression(italic("pha-4") * " enrichment, replicate 1"),
      ylab=expression(italic("pha-4") * " enrichment, replicate 2"))
  # XXX drawing axes separately, since we need font=5 for pch=183
  axis(1)
  axis(2)

  par(new=TRUE)
  plot(r[rownames(pha4.late), s1], r[rownames(pha4.late), s2],
    main="", xlab="", ylab="", col="#4040ffe0",
    pch=183, font=5, xlim=lim, ylim=lim, xaxt="n", yaxt="n")
  abline(0,1, col="#00000060")

  par(new=TRUE)
  plot(r[rownames(pha4.early), s1], r[rownames(pha4.early), s2],
    main="", xlab="", ylab="", col="#ff4040e0",
    pch=183, font=5, xlim=lim, ylim=lim, xaxt="n", yaxt="n")


  legend("topleft", c("Early", "Late"),
    pch=20, col=c("#ff0000a0", "#0000ffa0"))
}

pdf("git/sort_paper/plot/pha4KnownPharyngealScatterplot.pdf",
  width=4, height=4)   # was 5x5
par(mar=c(4,4,0.5,0.5)+0.1)
plot.scatter(pha4.replicates[1], pha4.replicates[2])
# plot.scatter(pha4.replicates[1], pha4.replicates[3])
# plot.scatter(pha4.replicates[2], pha4.replicates[3])

# p values for enrichments of known genes
r.pha4 = apply(r[,pha4.replicates], 1, mean)


if (FALSE) {

hist(r[ rownames(pha4.early) , "pha-4 9/1" ],
  main="Early pharyngeal genes", cex.main=0.8,
  xlab="enrichment", xlim=c(-2,6), col="grey")
label.panel("a)")
hist(r[ rownames(pha4.late) , "pha-4 9/1" ],
  main="Late pharyngeal genes", cex.main=0.8,
  xlab="enrichment", xlim=c(-2,6), col="grey")
label.panel("b)")

plot.it(r[ rownames(pha4.early), pha4.replicates ])
mtext("Early pharyngeal genes", cex=1.2, font=2)
legend("topleft", pha4.replicates, fill=replicate.colors,
  cex=0.8)
label.panel("a)", gp=gpar(fontsize=12, col="black"))
plot.it(r[ rownames(pha4.late), pha4.replicates ])
mtext("Late pharyngeal genes", cex=1.2, font=2)
label.panel("b)", gp=gpar(fontsize=12, col="black"))
}
dev.off()

#cat("early genes > 0 p =",
#  2 * t.test(r.pha4[rownames(pha4.early)], alternative="greater")$p.value)
#cat("late genes > 0 p =",
#  2 * t.test(r.pha4[rownames(pha4.late)], alternative="greater")$p.value)
cat("combined genes > 0 p =",
    t.test(r.pha4[c(rownames(pha4.early), rownames(pha4.late))], alternative="greater")$p.value, "\n")


