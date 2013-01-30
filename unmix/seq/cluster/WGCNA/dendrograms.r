# Plots dendrograms from WGCNA analysis.
# FIXME: combine these into one picture somehow?

load("git/unmix/seq/cluster/WGCNA/wnet/wnet.Rdata")

moduleColors = labels2colors(wnet$colors)

pdf("git/unmix/seq/cluster/WGCNA/dendrograms.pdf",
  title="WGCNA dendrograms", width=10, height=7.5)

for(i in 1:length(wnet$dendrograms))
plotDendroAndColors(wnet$dendrograms[[i]],
  moduleColors[wnet$blockGenes[[i]] ],
  "", main=paste("Gene dendrogram, block", i),
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05)

dev.off()

