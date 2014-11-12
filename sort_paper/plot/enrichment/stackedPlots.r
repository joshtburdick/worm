# Large "mashup" plots of what's enriched in
# FACS-sorted cells or clusters.

source("git/utils.r")

source("git/sort_paper/plot/enrichment/heatmapUtils.r")

source("git/sort_paper/tf/motifInfo.r")
source("git/tf/motif/motifName.r")

anatomy.info.matrix = function(f) {
  max.color.p = 6

  r = read.table(file=paste0("git/sort_paper/enrichment/summary/anatomyEnrichment/", f, ".tsv"),
    sep="\t", quote="", header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE)
  # XXX
  r = r[ r$group.name != "Tissue" & !grepl("depleted", r$cluster) , ]

  r = r[ -log10(r$p.corr) >= 1 , ]
   # & r$term.depth <= 4 , ]

  a = as.matrix(make.sparse.matrix(r$cluster, r$group.name, -log10(r$p.corr))) / max.color.p
  a = a[ , pick.top.few.columns(a, 1) ]
  a
}

gene.ontology.matrix = function(f) {
  p.cutoff = 1
  max.color.p = 6

  r = read.table(file=paste0("git/sort_paper/enrichment/geneOntology/", f, ".tsv"),
    sep="\t", quote="", header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE)
  r = r[ -log10(r$p.corr) >= p.cutoff & r$term.depth <= 4 
    & !grepl("depleted", r$cluster), ]
  a = as.matrix(make.sparse.matrix(r$cluster, r$Term, -log10(r$p.corr))) / max.color.p
  a = a[ , pick.top.few.columns(a, 2) ]
  a
}

motif.chip.matrix = function(f, p.cutoff=1, max.color.p=20, num.to.include=2) {
  r = read.tsv(paste0("git/sort_paper/tf/summary/",
    f, ".tsv"))

  r = r[ -log10(r$p.corr) >= p.cutoff & !grepl("depleted", r$group) , ]
  a = as.matrix(make.sparse.matrix(r$group, r$motif, -log10(r$p.corr))) / max.color.p
#  a[ a > 1e10 ] = 1e10
  a = a[ , pick.top.few.columns(a, num.to.include) ]

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
#   cluster.subset - list of clusters to include (optional)
# Side effects: plots a heatmap.
plot.stacked = function(f, cluster.subset = NULL) {

  anatomy.m = anatomy.info.matrix(f)
  go.m = gene.ontology.matrix(f)
  motif.m = motif.chip.matrix(paste0("motif/", f),
    p.cutoff=4, max.color.p = 7, num.to.include=1)
  chip.m = motif.chip.matrix(paste0("chip/", f),
    p.cutoff=1, max.color.p = 4, num.to.include=2)

  # XXX combine the rows appropriately
  cl = unique(c(rownames(anatomy.m), rownames(go.m),
    rownames(motif.m), rownames(chip.m)))

  # sort columns of these by clustering them
  sort.columns = function(a) {
    a[a > 1] = 1
    a[ , hclust(cor.dist(t(a)))$order ]
  }
  anatomy.m = sort.columns(anatomy.m)
  go.m = sort.columns(go.m)
  motif.m = sort.columns(motif.m)
  chip.m = sort.columns(chip.m)

  anatomy.m = add.rows(anatomy.m, cl)
  go.m = add.rows(go.m, cl)
  motif.m = add.rows(motif.m, cl)
  chip.m = add.rows(chip.m, cl)

  # compute "row" ordering
  r = cbind(anatomy.m, cbind(go.m, cbind(motif.m, chip.m)))

# browser()

  r[is.na(r)] = 0
  r[r > 1] = 1      # XXX deal with Infs
  cl = cl[ hclust(cor.dist(r))$order ]

  r = r[cl,]

  # possibly subset this
  if (!is.null(cluster.subset)) {
    r = r[ rownames(r) %in% cluster.subset , ]
  }

  # only keep cases in which something is present
  r = r[ , apply(r>0, 2, any) ]

#  color.scale = hsv(2/3, 0:255/255, 0.75)
  color.scale = hsv(0, 0, 255:0/255)

  r[ r==0 ] = NA
  rownames(r) = sub(" enriched", "", rownames(r))

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

    axis(2, at=(0:(dim(r)[2]-1)) / (dim(r)[2]-1), labels=ortho, las=2, cex.axis=0.3, line=1, tick=FALSE)

  # color different portions of the graph
#  rect(0, 0, 1, 1, border=NA,
#    col=hsv(0, 0.8, 1, alpha=0.2))
  color.columns = function(m, hue) {
    colnames.to.color = colnames(m)
    y = range(which(colnames(r) %in% colnames.to.color)) - 1
    n1 = length(colnames(r))-1
    rect(-0.5 / nrow(r), (y[1]-0.5)/n1,
      (nrow(r)+0.5)/nrow(r), (y[2]+0.5)/n1,
      border=hsv(0,0,0.6), lwd=0.2,      
      col=hsv(hue, 0.8, 1, alpha=0.15), xpd=TRUE)
  }

  color.columns(anatomy.m, 0)
  color.columns(go.m, 0.2)
  color.columns(motif.m, 0.4)
  color.columns(chip.m, 0.6)

  # ??? return info needed by highlight.column() ?
}

# Highlights one column of the graph.
highlight.column = function(colnames, a, hue) {
  x1 = which(colnames == a)[1]
  x = c(x1-1.5, x1-0.5) / (length(colnames)-1)
# browser()
  rect(x[1], -5, x[2], 5, lwd=0.3,
    border=hsv(hue, 1, 0.6, alpha=0.5), col=hsv(hue, 0.8, 1, alpha=0.2))
}

system(paste("mkdir -p git/sort_paper/plot/enrichment/stackedPlots"))

# a subset of the clustering
pdf("git/sort_paper/plot/enrichment/stackedPlots/hier.300.subset1.pdf",
  width=5, height=10)
par(mar=c(6,15,1,1))
ao = anatomy.info.matrix("hier.300.clusters")
plot.stacked("hier.300.clusters", rownames(ao))
     # , as.character(c(1,2,30,52,79,223,286)))
dev.off()

pdf("git/sort_paper/plot/enrichment/stackedPlots/hier.300.pdf",
  width=16.5, height=14.5)
par(mar=c(6,15,1,1))

plot.stacked("hier.300.clusters")

# possibly interesting clusters
# highlight.column(rownames(r), "52", 0)
# highlight.column(rownames(r), "110", 1/4)
# highlight.column(rownames(r), "286", 2/4)

dev.off()



# things enriched in FACS-sorted fractions
pdf("git/sort_paper/plot/enrichment/stackedPlots/facs.pdf",
  width=4.5, height=6.8)
par(mar=c(6,15,1,1))
plot.stacked("facs")
dev.off()

