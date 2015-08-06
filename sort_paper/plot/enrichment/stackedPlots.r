# Large "mashup" plots of what's enriched in
# FACS-sorted cells or clusters.

source("git/utils.r")
source("git/plot/utils.r")
source("git/sort_paper/plot/enrichment/heatmapUtils.r")

# source("git/tf/motif/motifName.r")
source("git/sort_paper/tf/motif/hughes/motifInfo.r")
source("git/sort_paper/tf/motifInfo.r")
source("git/data/wormbase/wb.cluster.name.r")

# clusters to annotate in the "selected clusters" thing
selected.clusters = as.matrix(read.tsv(
  "git/sort_paper/plot/enrichment/selectedClusters.tsv"))[,1]

# Reads anatomy info (or WB cluster info).
anatomy.info.matrix = function(f, ao.or.wbcluster) {
  max.color.p = 9

  r = read.table(file=paste0("git/sort_paper/enrichment/summary/",
    ao.or.wbcluster, "/", f, ".tsv"),
    sep="\t", quote="", header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE)
  # XXX
  r = r[ r$group.name != "Tissue" , ]

  r = r[ -log10(r$p.corr) >= 1 , ]    # was 1
   # & r$term.depth <= 4 , ]

  group.names = if (ao.or.wbcluster=="anatomyEnrichment")
    r$group.name
  else
    r$group

  a = as.matrix(make.sparse.matrix(r$cluster, group.names, -log10(r$p.corr))) / max.color.p
  a = a[ , pick.top.few.columns(a, 1) ]
  a
}

gene.ontology.matrix = function(f) {
  p.cutoff = 1
  max.color.p = 8

  r = read.table(file=paste0("git/sort_paper/enrichment/geneOntology/", f, ".tsv"),
    sep="\t", quote="", header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE)
  r = r[ -log10(r$p.corr) >= p.cutoff & r$term.depth <= 4 , ]
  a = as.matrix(make.sparse.matrix(r$cluster, r$Term, -log10(r$p.corr))) / max.color.p
  a = a[ , pick.top.few.columns(a, 2) ]
  a
}

motif.chip.matrix = function(f, p.cutoff=1, max.color.p=20, num.to.include=2) {
  r = read.tsv(gzfile(paste0("git/sort_paper/tf/motif/hyperg/",
    f, ".tsv.gz")))
  colnames(r)[[1]] = "motif"

  r = r[ -log10(r$p.corr) >= p.cutoff , ]
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

# Draws a p-value scale.
draw.scale = function(dims, x) function(hue, y, max.p) {
  d = dims - 1

#  image(c(x-1,x) / d[1],
#    (y:(y+4)) / d[2],
#    z = t(as.matrix(0:4/4)), add=TRUE, useRaster=FALSE,
#    xaxs="n", yaxs="n", xpd=FALSE)
  height = 7
  n = 10
  rect((x-1) / d[1], (y + c(0:(n-1)) * (height/n)) / d[2],
    x / d[1], (y + c(1:n) * (height/n)) / d[2],
    col = hsv(0,0,c(n:1)/n), lty=0, xpd=NA)
  rect((x-1) / d[1], y / d[2],
    x / d[1], (y + height) / d[2],
    col = hsv(hue, 0.8, 1, alpha=0.15), lty=0, xpd=NA)

  text(c(x+1, x+1) / d[1], c(y, y+height) / d[2],
    labels=c(0, max.p), xpd=NA, cex=0.4)
  text((x-2.5) / d[1], y / d[2],
    labels=expression("  -log"[10] * " p"),
    xpd=NA, cex=0.5, srt=90, adj=c(0,1))
}

# For each motif, a concise list of potential orthologs
# (prettyprinted).
ortholog.by.motif.prettyprint = {

  # Possibly italicizes a gene name, as an expression.
  format.gene = function(g) {
    if (grepl("^[a-z][a-z][a-z][a-z]?-[0-9][0-9]?[0-9]?",
        ignore.case=FALSE, g))
      italicize(g)
    else
      as.expression(g)
  }

  # first, the name of whatever the motif was
  a = list()
  org.name.short = c("C.elegans" = "Ce", "D.melanogaster" = "Dm",
    "H.sapiens" = "Hs", "M.musculus" = "Mm")
  for(m in motif.info$motif.id) {

    org = org.name.short[ motif.info[m, "species" ] ]
    related.gene = motif.info[ m, "related.gene" ]

    # italicize gene name if it's Dm, or a known C.elegans gene
    if (motif.info[m,"species"] == "D.melanogaster" ||
      (motif.info[m,"species"] == "C.elegans" &&
        grepl("^[a-z][a-z][a-z][a-z]?-[0-9][0-9]?[0-9]?",
          ignore.case=FALSE, related.gene))) {
      a[[m]] = expr.format(expression(bold(species) * " " * bolditalic(gene)),
        list(species = org, gene = related.gene))
    }
    else {
      a[[m]] = expr.format(expression(bold(species) * " " * bold(gene)),
        list(species = org, gene = related.gene))
    }

if (FALSE) {
    x = expression(1 * 2)
    x[[1]][[2]] = (paste(org.name.short[ motif.info[ m, "species" ] ], " "))
    x[[1]][[3]] = (motif.info[ m, "related.gene" ])
    if (motif.info[m,"species"] == "D.melanogaster")
      x[[1]][[3]] = italicize(motif.info[m,"related.gene"])[[1]]
    if (motif.info[m,"species"] == "C.elegans" &&
        grepl("^[a-z][a-z][a-z][a-z]?-[0-9][0-9]?[0-9]?",
          ignore.case=FALSE, motif.info[m,"related.gene"]))
      x[[1]][[3]] = italicize(motif.info[m,"related.gene"])[[1]]
    a[[m]] = x
}

  }

  # Concatenates gene names
  expr.concat = function(s) {
    x = format.gene(s[1])

    if (length(s) > 1) {
      for(i in 2:length(s)) {
        x1 = expression(1*" "*2)
        x1[[1]][[2]][[2]] = x[[1]]
        x1[[1]][[3]] = format.gene(s[i])[[1]]
        x = x1
      }
    }
    return(x)
  }

  for(m in names(a)) {
    g = unique(as.character(motif.ortholog[ motif.ortholog$motif.id==m, "gene" ]))

    if (length(g) > 0) {
      # slight reordering
      g = g[ order( !(g %in% c("pha-4"))) ]

      # if there are many genes, summarize this
      if (length(g) > 5) {
        g = c(g[1:4], paste("and", length(g)-4, "others"))
      }

      x1 = expression(1 * " (" * 2 * ")")
      x1[[1]][[2]][[2]][[2]] = a[[m]][[1]]
      x1[[1]][[2]][[3]] = expr.concat(g)[[1]] 
      a[[m]] = x1    #    paste(g, collapse = " ")   
    }
  }

  a   # c(a, recursive=TRUE)
}

# Plots one of these heatmaps.
# Args:
#   f - name of file
#   cluster.subset - vector with names = clusters to include
#     (optional). (Entries of this will be labels.)
# Side effects: plots a heatmap.
plot.stacked = function(f, cluster.subset = NULL) {

  # read various enrichments
  anatomy.m = anatomy.info.matrix(f, "anatomyEnrichment")
  cluster.m = anatomy.info.matrix(f, "wormbaseCluster")
  go.m = gene.ontology.matrix(f)

  # read motif / ChIP signal enrichment
  # XXX filenames not yet changed over
  f.chip = f
  f.motif = f
  if (f == "facs") {
    f.chipb = "facs_enrichmentVsOpposite_2"
    f.motif = "facs_vs_opposite_1"
  }

  motif.m = motif.chip.matrix(paste0("summary/hughes/", f.motif),
    p.cutoff=0, max.color.p = 10, num.to.include=1)
  chip.m = motif.chip.matrix(paste0("chipTable/", f.chip),
    p.cutoff=0, max.color.p = 5, num.to.include=1)

  # filter rows of m to be non-redundant by motif gene name
  j = ! duplicated(motif.gene[colnames(motif.m)])
  motif.m = motif.m[ , j ]
  cat("reduced from", length(j), "to", sum(j), "rows\n")

  # possibly subset these
  if (!is.null(cluster.subset)) {
    f = function(a, num.to.keep) {
      a = row.subset(a, names(cluster.subset))
      a = a[ , pick.top.few.columns(a, num.to.keep) ]
      a
    }

    anatomy.m = f(anatomy.m, 1)
    cluster.m = f(cluster.m, 1)
    go.m = f(go.m, 2)
    motif.m = f(motif.m, 1)
    chip.m = f(chip.m, 1)
  }

  # sort columns of these by clustering them
  sort.columns = function(a) {
    if (ncol(a) < 2) {
      return(a)
    }
    a[a > 1] = 1
    a[ , hclust(cor.dist(t(a)))$order ]
  }
  anatomy.m = sort.columns(anatomy.m)
  cluster.m = sort.columns(cluster.m)
  go.m = sort.columns(go.m)
  motif.m = sort.columns(motif.m)
  chip.m = sort.columns(chip.m)

  # some renaming
  colnames(cluster.m) = sapply(colnames(cluster.m), cluster.name.format)
  colnames(chip.m) = gsub("_", " ", colnames(chip.m))

  # XXX combine the rows appropriately
  cl = unique(c(rownames(anatomy.m), rownames(cluster.m),
    rownames(go.m), rownames(motif.m), rownames(chip.m)))

  anatomy.m = add.rows(anatomy.m, cl)
  cluster.m = add.rows(cluster.m, cl)
  go.m = add.rows(go.m, cl)
  motif.m = add.rows(motif.m, cl)
  chip.m = add.rows(chip.m, cl)

  # compute "row" ordering
#  r = cbind(anatomy.m, cbind(cluster.m, cbind(go.m, cbind(motif.m, chip.m))))

  # XXX hack
  e = 1e-100
  r = cbind(anatomy.m, e, cluster.m, e, go.m, e, motif.m, e, chip.m)
  colnames(r) = gsub("^e$", "", colnames(r))

  r[is.na(r)] = 0
  r[r > 1] = 1      # XXX deal with Infs
  cl = cl[ hclust(cor.dist(r))$order ]

  r = r[cl,]

  # only keep cases in which something is present
  r = r[ , apply(r>0, 2, any) ]

  # XXX hack to show ordering by "number of genes in group"
  if (FALSE) {
    source("git/sort_paper/FACS/enrichedInFraction.r")
    sort.by.count = sort(sapply(facs.enriched.depleted, sum))
    sort.by.count = sort.by.count[ sort.by.count > 0 ]
    r = r[ names(sort.by.count) , ]
  }

#  color.scale = hsv(2/3, 0:255/255, 0.75)
  color.scale = hsv(0, 0, 255:0/255)

  r[ r==0 ] = NA
#  rownames(r) = sub(" enriched", "", rownames(r))

  # possibly add to cluster subset names
  if (!is.null(cluster.subset)) {
    r1 = cluster.subset[as.character(rownames(r))]
    rownames(r) = ifelse(r1 == "", rownames(r),
      paste0(rownames(r), ": ", r1))
  }

  image(r, col=color.scale, xaxt="n", yaxt="n", bty="n", zlim=c(0,1))
  axis(1, at=(0:(dim(r)[1]-1)) / (dim(r)[1]-1), labels=rownames(r), las=2, cex.axis=0.3, line=-0.9, tick=FALSE)

  rownames1 = colnames(r)

if (FALSE) {
  # orthology info
  ortho = rep("", length(rownames1))
  names(ortho) = rownames1
  m = rownames1 %in% names(ortholog.by.motif.small)
  ortho[ rownames1[ m ] ] = ortholog.by.motif.small[ rownames1 [ m ] ]

  # convert motif names to the relevant gene name
  m = rownames1 %in% names(motif.gene)
  rownames1[ m ] = motif.gene[ rownames1[m] ]
}

  # combined motif names and ortholog info
  m = rownames1 %in% names(ortholog.by.motif.prettyprint)
  rownames1[ m ] = ortholog.by.motif.prettyprint[ rownames1[ m ] ]
# browser()

  # label rows
  axis(2, at=(0:(dim(r)[2]-1)) / (dim(r)[2]-1),
    labels=as.expression(c(rownames1, recursive=TRUE)),
    las=2, cex.axis=0.3, line=-0.9, tick=FALSE)

  # label orthologs (currently disabled)
  if (FALSE) {
    axis(2, at=(0:(dim(r)[2]-1)) / (dim(r)[2]-1),
      labels=sapply(ortho, italicize),
      las=2, cex.axis=0.3, line=1.5, tick=FALSE)
  }

  # add grid lines
  g = dim(r) - 1
  abline(h = (c(0:g[2])-0.5) / g[2], lwd=0.1, col="#00000040")
  # vertical lines will be added later

  # color different portions of the graph, and label them
#  rect(0, 0, 1, 1, border=NA,
#    col=hsv(0, 0.8, 1, alpha=0.2))
  color.columns = function(m, hue, ylab, max.p) {
    colnames.to.color = colnames(m)
    y = range(which(colnames(r) %in% colnames.to.color)) - 1
    n1 = length(colnames(r))-1
    rect(-0.5 / nrow(r), (y[1]-0.5)/n1,
      (nrow(r)+0.5)/nrow(r), (y[2]+0.5)/n1,
      border=hsv(0,0,0.6), lwd=0.2,      
      col=hsv(hue, 0.8, 1, alpha=0.15), xpd=TRUE)

    for(x in (c(0:g[1])-0.5) / g[1]) {
      lines(c(x, x), c(y[1]-0.5, y[2]+0.5)/n1, lwd=0.1, col="#00000040")
    }

    text(-30/nrow(r), 0.5*(y[1]+y[2]) / ncol(r),
      labels=ylab,
      srt=90, xpd=TRUE, adj=0.5, cex=0.68)

    # XXX hack
    if (FALSE && ylab=="Possible orthologs") {
      text(-7/nrow(r), 0.5*(y[1]+y[2]) / ncol(r), 
        labels="Motifs",
        srt=90, xpd=TRUE, adj=0.5, cex=0.68)
    }
  }

  color.columns(anatomy.m, 0, "Anatomy\nterms", 9)
  color.columns(cluster.m, 0.2, "Expression\nclusters", 9)
  color.columns(go.m, 0.4, "GO terms\n", 8)
  color.columns(motif.m, 0.6, "Motifs\n", 10)
  color.columns(chip.m, 0.8, "ChIP\nsignals", 5)

  # ??? return info needed by highlight.column() ?

  # add on a p-value scale
  draw.scale(dim(r), -23.5)(0, -10, 10)
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

if (TRUE) {

# a subset of the clustering
pdf("git/sort_paper/plot/enrichment/stackedPlots/hier.300.subset1.pdf",
  width=3, height=6)
par(mar=c(5,10,0.1,0.1))
ao = anatomy.info.matrix("hier.300.clusters", "anatomyEnrichment")
wbc = anatomy.info.matrix("hier.300.clusters", "wormbaseCluster")
# cl.subset = unique(c(rownames(ao), sort(rownames(wbc), decreasing=TRUE)[1:5]))

plot.stacked("hier.300.clusters", selected.clusters)
     # , as.character(c(1,2,30,52,79,223,286)))
mtext("Cluster", side=1, line=3.5, cex=0.68)

dev.off()

# all of the clustering
pdf("git/sort_paper/plot/enrichment/stackedPlots/hier.300.pdf",
  width=21.5, height=24)
par(mar=c(3,11,0.1,0.1))

plot.stacked("hier.300.clusters")

mtext("Cluster", side=1, line=1.5, cex=0.68)
dev.off()
}

# things enriched in FACS-sorted fractions
pdf("git/sort_paper/plot/enrichment/stackedPlots/facs.pdf",
  width=4.2, height=6.5)
par(mar=c(5,10,0.1,0.1))
plot.stacked("facs")
mtext("Sort fraction", side=1, line=4, cex=0.68)
dev.off()

