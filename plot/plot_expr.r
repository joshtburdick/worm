# Utilities for drawing expression as a tree.

# source("tree/two_color_image.r")        # may be omitted
# source("tree/deconvolve_tree.r")
# source("/murrlab/jburdick/gcb/work/svnroot/trunk/deconvolve/R/tree/lineage_tree.r")

source("R/lineage/tree_utils.r")

cell.time.on.off =
  read.csv(gzfile("data/image/cellTimeOnOff.csv.gz"),
    row.names=1)
# XXX fix EMS mis-ordering
cell.time.on.off = cell.time.on.off[ lin.node.names , ]

# Default list of lineages to label.
get.interior.nodes = function() {
  int.n = lin.node.names[ nchar(lin.node.names) <= 3 ]
  int.n = union(int.n, grep("AB...$|", lin.node.names, value=TRUE))

  int.n = grep("AB.?.?.?.?$|MS.?.?$|E.?$|C.?$|D.?$|P1|P2|P3", lin.node.names, value=TRUE)

  int.n
}
int.n = get.interior.nodes()

# Converts an expression dataset into a
# series of line segments, drawable using segments().
make.cell.time.segments = function(cell.time) {

  ct1 = cell.time
  
  x0 = cell.to.column[as.vector(ct1$cell)]
  y0 = as.numeric(ct1$time)
  x1 = x0
  y1 = y0 + 1
  
  list(x0=x0, y0=y0, x1=x1, y1=y1)
}

# Makes just the horizontal lines from the first node at a time
# to its parent.
make.horiz.segments = function(cell.time) {

  ct1 = cell.time

  n = as.vector(names(parent.of)) 
  x0 = cell.to.column[ as.vector(parent.of) ]
  y0 = onset.time.1[ as.vector(parent.of) ]
  x1 = cell.to.column[ as.vector(names(parent.of)) ]
  y1 = y0
  
# currently doesn't add color.
#  h.col = paste(as.vector(parent.of),
#    onset.time(as.vector(parent.of)), sep=":")

  list(x0=x0, y0=y0, x1=x1, y1=y1)
}

# Scales positive data to [0,1].
# (Negative numbers are clipped to 0, and NAs are set to 0.)
# FIXME use something more robust than "min" and "max" ?
scale.to.unit = function(x) {
  x[ is.na(x) ] = 0
  x[ x < 0 ] = 0
  x1 = x - min(x, na.rm=TRUE)
  x2 = x1 / max(x1, na.rm=TRUE)

  # XXX presumably these aren't needed...
  x2[ x2 < 0 ] = 0
  x2[ x2 > 1 ] = 1

  x2[ is.na(x2) ] = 0   # this may not matter
  x2
}

# Given a dendrogram, computes columns such that the leaves
# are all evenly spaced, and intermediate nodes are halfway
# in between those.
number.lineage.columns = function(lin) {

  # XXX using this instead of "c(..., recursive=TRUE)",
  # which isn't working, for what reason I know not
  get.leaf.nodes = function(x) {
    if (attr(x, "leaf"))
      return(attr(x, "label"))

    return(c(get.leaf.nodes(x[[1]]),
      get.leaf.nodes(x[[2]])))
  }

  # start with where leaves are located
  node.names = names(lineage.to.list(lin))
  node.to.index = rep(NA, length(node.names))
  names(node.to.index) = node.names
  leaves = get.leaf.nodes(lin)       # was c(lin, recursive=TRUE)
  node.to.index[leaves] = 1:length(leaves)

  # number internal nodes, by traversing the tree in post-order
  walk.tree = function(l) {
    name = attr(l, "label")
    if( !attr(l, "leaf")) {
      walk.tree(l[[1]])
      walk.tree(l[[2]])
      node.to.index[name] <<-
        (node.to.index[ attr(l[[1]], "label") ] +
        node.to.index[ attr(l[[2]], "label") ]) / 2
    }
  }
  walk.tree(lin)

  node.to.index
}

# Plots an expression data set, including time.
# Args:
#   r - data frame with columns
#     cell - cell name (as a string, not a factor)
#     time.1, time.2 - start and end time for this segment (in minutes)
#     col - a color (as a string; see rgb() or hsv())
#   main - title
#   root (= P0) - root of the tree
#   times (= c(0, 350) - time interval to include
#   lwd - width of lines
# XXX assumes Sulston-style node naming.
# Side effects: draws a tree, colored according to the given colors.
plot.segments = function(r, main, root="P0",
  times=c(0, 350), lwd=3) {

  # part of the lineage to plot
  lin1 = get.lineage.by.name(lin, root)

  num.leaves = length(c(lin1, recursive=TRUE))
  r = r[ r$cell %in% names(lineage.to.list(lin1)) , ]

  # compute x coordinates, based on leaf positions
  cell.to.x = number.lineage.columns(lin1)
  r$x = cell.to.x[r$cell]

  # get first and last times for each cell
  cells = unique(r$cell)

  time.1 = c(by(r$time.1, r$cell, min))
  time.2 = c(by(r$time.2, r$cell, max))

  xlim=c(0,num.leaves)
  ylim=c(times[2], times[1])

  # set up graph
  par(bg="#ffffff")
  plot(0, 0, xlim=xlim, ylim=ylim,
    main=main, xlab="", ylab="time", xaxt="n", type="n")
print(xlim)

  # plot the tree: first one line per cell, then the connecting line
  segments(cell.to.x[cells], time.1[cells], cell.to.x[cells], time.2[cells],
    col="#c0c0c0", xlim=xlim, ylim=ylim, lwd=lwd, lend=2)
  a = lineage.triples(lin)
  segments(cell.to.x[a$daughter.1], time.1[a$daughter.1],
    cell.to.x[a$daughter.2], time.1[a$daughter.2],
    col="#c0c0c0", xlim=xlim, ylim=ylim, lwd=lwd, lend=2)

  # then, plot the data
  segments(cell.to.x[r$cell], r$time.1, cell.to.x[r$cell], r$time.2,
    col=r$col, xlim=xlim, ylim=ylim, lwd=lwd, lend=2)

  # add labels for a few branches in the tree
  int.n.1 = int.n  # %in% lineage.to.list(lin1)
  # XXX hack
  if (root == "ABplp") {
    int.n.1 = c(int.n.1, "ABplpaa", "ABplpap", "ABplppa", "ABplppp")
  }
  text(cell.to.x[int.n.1], time.1[int.n.1]-5, int.n.1,
    adj=0, cex=0.7, srt=90)

#  axis(1, at=cell.to.column[axis.nodes], labels=axis.nodes, cex=0.5, las=3)

}

# Plots expression data sets as red, green, and blue.
# Assumes that the expression data vectors' lengths are the
# the same as the number of rows in cell.time .
# Possibly not used...
plot.segments.old = function(cell.time, main, x.red, x.green, x.blue) {

  # set up coordinates
  par(bg="#ffffff")
  plot(0,0,xlim=c(0,length(lin.node.names)), ylim=c(551,0),
    main=main, xlab="", ylab="time", xaxt="n", type="n")

  # plot the tree. XXX this is slow
  seg.1 = make.cell.time.segments(cell.time.all)
  segments(seg.1$x0, seg.1$y0, seg.1$x1, seg.1$y1, col="#707070", lwd=1)
  seg.2 = make.horiz.segments(cell.time.all)
  segments(seg.2$x0, seg.2$y0, seg.2$x1, seg.2$y1, col="#707070", lwd=1)

  # then, plot the data
  seg = make.cell.time.segments(cell.time)
  col=rgb(scale.to.unit(x.red), scale.to.unit(x.green), scale.to.unit(x.blue))
  segments(seg$x0, seg$y0, seg$x1, seg$y1, col=col, lwd=5)
  
  # add labels for a few branches in the tree
  int.n = lin.node.names[ nchar(lin.node.names) <= 3 ]
  text(cell.to.column[int.n], onset.time[int.n]-10, int.n, cex=0.9)
  axis.nodes = lin.node.names[ nchar(lin.node.names) <= 5 ]
  axis(1, at=cell.to.column[axis.nodes], labels=axis.nodes, cex=0.5, las=3)
}

# Plots per-cell expression data.
# Args:
#   col - the colors, as a vector with names the same as the
#     "standard" 1341 cells.
#   cell.time.on.off - matrix indicating which cell is on when
plot.segments.per.cell = function(col, cell.time.on.off,
    times=c(-30, 350), lwd=5) {

  ylim = c(times[2], times[1])

  # set up coordinates
  par(bg="#ffffff")
  plot(0,0,xlim=c(1,length(lin.node.names)), ylim = ylim,  # ylim=c(551,0),
    main="", xlab="", ylab="time (minutes)", xaxt="n", type="n")
  cell.x = cell.to.column[ rownames(cell.time.on.off) ]

  # plot the tree
  segments(cell.x, cell.time.on.off$on, cell.x, cell.time.on.off$off,
    col="#707070", lwd=lwd, lend=2)
  # FIXME draw the tops of trees here
  segments(cell.to.column[ as.vector(rownames(cell.time.on.off)) ],
    cell.time.on.off$on,
    cell.to.column[ parent.of[ as.vector(rownames(cell.time.on.off)) ]],
    cell.time.on.off$on,
    col="#d0d0d0", lwd=lwd, lend=2)   # was #707070, lwd=1

  # then, plot the data
  segments(cell.x, cell.time.on.off$on, cell.x, cell.time.on.off$off,
    col=col[rownames(cell.time.on.off)], lwd=lwd, lend=2)

  # add labels. FIXME label this better?
#  axis.nodes = lin.node.names[ nchar(lin.node.names) <= 5 ]
#  axis(1, at=cell.to.column[axis.nodes], labels=axis.nodes, cex=0.5, las=3)

  # add labels for a few branches in the tree
#  int.n = lin.node.names[ nchar(lin.node.names) <= 3 ]
#  text(cell.to.column[int.n], cell.time.on.off[int.n,"on"]-10,
#    int.n, cex=0.9)
  text(cell.to.column[int.n], cell.time.on.off[int.n,"on"]-5,
    int.n, cex=0.7, adj=0, srt=90)
}

# Plots where different groups of reporters share the same value.
# XXX moving this elsewhere
# Args:
#   m - a matrix with one row per reporter
#   num.reporters - number of reporters to include
# Side effects: writes a PNG showing where there are different
#   reporter combinations.
# XXX broken: colors seem to be weirdly quantized
# also deprecated, as I'm using per-cell info (for now)
plot.reporter.locations = function(m, reporters) {
  num.reporters = length(reporters)
  m1 = t( m[reporters,] )
  gr = group.m(m1)

  name = paste(num.reporters, "reporters")
  png(paste(name, ".png", sep=""), width=3000, height=600)

#  colors = hsv(runif(max(gr$row.map)), 0.9, 0.3)
  colors = hsv(c(1:max(gr$row.map)) /max(gr$row.map), 0.9, 0.9)
  col = make.numbering(colors)[ gr$row.map ]

  # FIXME XXX segments() should probably expect the list of colors
  col.rgb = col2rgb(col)
cat("length(col) =", length(col), "\n")
cat("length(col.rgb[red,]) =", length(col.rgb["red",]), "\n")
  r1 = col.rgb["red",]
  g1 = col.rgb["green",]
  b1 = col.rgb["blue",]
  names(r1) = rownames(cell.time)
  names(g1) = rownames(cell.time)
  names(b1) = rownames(cell.time)

  plot.segments(cell.time, name, r1, g1, b1)

  dev.off()
}





###### the following are various examples of calling plot.segments.per.cell,
# and are deprecated

# Plot actual expression patterns; this is largely to test the
# tree-drawing code.
# Args:
#   expr - the expression matrix
#   output.subdir - subdirectory (in "out") in which to write images
plot.expr.all = function(expr, output.subdir) {
  system(paste("mkdir -p out/", output.subdir, sep=""))
  for(i in 1:dim(expr)[1]) {
#   for(i in 1:5) {
    name = rownames(expr)[i]
cat(name, "")
    expr1 = expr[i,]
#    expr.is.na = is.na(expr1)
#    expr1[ is.na(expr1) ] = 0
    expr1[ expr1 <= 0 ] = 0

    png(file=paste("out/", output.subdir, "/", name, ".png", sep=""),
      width=3000, height=600)

#    plot.segments(cell.time.all, name, expr1, 0*expr1, expr.is.na)

    col = rgb(scale.to.unit(expr1), 0, 0)
    names(col) = names(expr1)
    plot.segments.per.cell(col)
    title(name)
    dev.off()
  }

cat("\n")
}

# Plot actual expression patterns; this is largely to test the
# tree-drawing code.
# Args:
#   x.red, x.green, x.blue - the expression matrix
#   output.subdir - subdirectory (in "out") in which to write images
plot.expr.all.multicolor = function(x.red, x.green, x.blue, output.subdir) {
  system(paste("mkdir -p out/", output.subdir, sep=""))
  for(i in 1:dim(x.red)[1]) {
#   for(i in 1:5) {
    name = rownames(x.red)[i]
cat(name, "")

    png(file=paste("out/", output.subdir, "/", name, ".png", sep=""),
      width=3000, height=600)

    col = rgb(scale.to.unit(x.red[i,]),
      scale.to.unit(x.green[i,]),
      scale.to.unit(x.blue[i,]))
    names(col) = names(x.red[i,])
    plot.segments.per.cell(col)
    title(name)
    dev.off()
  }

cat("\n")
}

# Plot actual expression patterns, in a PDF file.
# Args:
#   x.red, x.green, x.blue - the expression matrix
#   output.file - PDF file in which to write images
plot.expr.all.multicolor.pdf = function(x.red, x.green, x.blue, output.file) {

  pdf(file=output.file, width=7.5, height=10)
  par(mfrow=c(3,1))

  for(i in 1:dim(x.red)[1]) {
    name = rownames(x.red)[i]
cat(name, "")

    col = rgb(scale.to.unit(x.red[i,]),
      scale.to.unit(x.green[i,]),
      scale.to.unit(x.blue[i,]))
    names(col) = names(x.red[i,])
    plot.segments.per.cell(col)
    title(name)
  }

  dev.off()
  cat("\n")
}

# Plots "actual" expression in red, and two different predictions in green;
# used for comparing methods.
# Args:
#   expr.ref, expr.1, expr.2 - expression datasets, with one column per cell.
#   output.subdir - subdirectory (in "out") in which to write images
plot.expr.comparison.red.green = function(expr.ref, expr.1, expr.2, output.subdir) {
  system(paste("mkdir -p out/", output.subdir, sep=""))
  for(i in 1:dim(expr.ref)[1]) {
#   for(i in 1:5) {
    name = rownames(expr.ref)[i]
cat(name, "")

    png(file=paste("out/", output.subdir, "/", name, ".png", sep=""),
      width=3000, height=1200)
    par(mfrow=c(2,1))

    # first dataset
    col = rgb(scale.to.unit(expr.ref[i,]),
      scale.to.unit(expr.1[i,]),
      scale.to.unit(0*expr.ref[i,]))
    names(col) = names(expr.ref[i,])
    plot.segments.per.cell(col)
    title(name)

    # second dataset
    col = rgb(scale.to.unit(expr.ref[i,]),
      scale.to.unit(expr.2[i,]),
      scale.to.unit(0*expr.ref[i,]))
    names(col) = names(expr.ref[i,])
    plot.segments.per.cell(col)
    title(name)

    dev.off()
  }

cat("\n")
}

