# Heatmap of eigengenes.

source("git/unmix/seq/cluster/WGCNA/eigengene.r")

# Draws a heatmap of non-zero coefficients
draw.heatmap = function(x, z.max=1) {

  red.blue.colors = c(hsv(2/3, 128:0/128, 1), hsv(0, 0:128/128, 1))

  if (is.null(z.max)) {
    z.max = max(abs(x), na.rm=TRUE)
  }

  x[ x > z.max ] = z.max
  x[ x < -z.max ] = -z.max

  image(x, xaxt="n", yaxt="n",
    col=red.blue.colors, zlim=c(-z.max, z.max),
    main="")
  axis(1, at=(0:(dim(x)[1]-1)) / (dim(x)[1]-1),
    labels=rownames(x), las=2, cex.axis=0.5)
  axis(2, at=(0:(dim(x)[2]-1)) / (dim(x)[2]-1),
    labels=colnames(x), las=2, cex.axis=0.5)
}

pdf("git/unmix/seq/cluster/WGCNA/eigengeneHeatmap.pdf",
  title="WGCNA eigengenes", width=7, height=7)

draw.heatmap(t(module.eigengene[c(64:1),]), z.max = 0.5)

dev.off()

