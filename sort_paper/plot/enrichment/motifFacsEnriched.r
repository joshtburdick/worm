# Large "mashup" plots of what's enriched in
# FACS-sorted cells or clusters.

source("git/utils.r")
source("git/sort_paper/plot/enrichment/heatmapUtils.r")
source("git/sort_paper/tf/motifInfo.r")
source("git/tf/motif/motifName.r")
source("git/data/wormbase/wb.cluster.name.r")

load("git/sort_paper/tf/motif/


motif.chip.matrix = function(f, p.cutoff=1, max.color.p=20, num.to.include=2) {
  r = read.tsv(gzfile(paste0("git/sort_paper/tf/motif/hyperg/",
    f, ".tsv.gz")))
  colnames(r)[[1]] = "motif"

  r = r[ -log10(r$p.corr) >= p.cutoff , ]
  a = as.matrix(make.sparse.matrix(r$group, r$motif, -log10(r$p.corr))) / max.color.p

# disabling this for now
#  a = a[ , pick.top.few.columns(a, num.to.include) ]

  a
}

# Hack to add blank rows to a matrix.
add.rows = function(a, rows) {
  b = matrix(nrow = length(rows), ncol = ncol(a))
  rownames(b) = rows
  colnames(b) = colnames(a)
  b[ rownames(a), colnames(a) ] = a
  b
}

# Plots one of these heatmaps.
# Args:
#   f - name of file
# Side effects: plots a heatmap.
plot.1 = function(f) {

  motif.m = motif.chip.matrix(paste0("table/", f),
    p.cutoff=0, max.color.p = 40)

  # sort columns of these by clustering them
  sort.columns = function(a) {
    if (ncol(a) < 2) {
      return(a)
    }
    a[a > 1] = 1
    a[ , hclust(cor.dist(t(a)))$order ]
  }

  motif.m = sort.columns(motif.m)

  # some renaming
#  colnames(motif.m) = sapply(colnames(motif.m), cluster.name.format)

  # compute "row" ordering
  r = motif.m

  r[is.na(r)] = 0
  r[r > 1] = 1      # XXX deal with Infs
  cl = rownames(r)
  cl = cl[ hclust(cor.dist(r))$order ]

  r = r[cl,]

  color.scale = hsv(0, 0, 255:0/255)

  r[ r==0 ] = NA
#  rownames(r) = sub(" enriched", "", rownames(r))
print(dim(r))
  image(r, col=color.scale, xaxt="n", yaxt="n", bty="n", zlim=c(0,1))
  axis(1, at=(0:(dim(r)[1]-1)) / (dim(r)[1]-1), labels=rownames(r), las=2, cex.axis=0.3, line=-0.9, tick=FALSE)

  rownames1 = colnames(r)

  # orthology info
  ortho = rep("", length(rownames1))
  names(ortho) = rownames1
  m = rownames1 %in% names(motif.info)
  ortho[ rownames1[ m ] ] = motif.info[ rownames1 [ m ] ]

  m = rownames1 %in% names(motif.gene)
  rownames1[ m ] = motif.gene[ rownames1[m] ]
  m = rownames1 %in% names(motif.name)
  rownames1[ m ] = motif.name[ rownames1[m] ]

    axis(2, at=(0:(dim(r)[2]-1)) / (dim(r)[2]-1), labels=rownames1, las=2, cex.axis=0.3, line=-0.9, tick=FALSE)

    axis(2, at=(0:(dim(r)[2]-1)) / (dim(r)[2]-1), labels=ortho, las=2, cex.axis=0.3, line=4, tick=FALSE)

}

# things enriched in FACS-sorted fractions
pdf("git/sort_paper/plot/enrichment/motifFacsEnriched.pdf",
  width=5.5, height=100)
par(mar=c(5,15,0.1,0.1))
plot.1("facs_0.892")
dev.off()

