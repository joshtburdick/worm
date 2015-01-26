# Writes out a clustered dataset in TreeView-readable
# format. Genes in each cluster will be hierarchically
# clustered together, and colored some color.

library(amap)
library(ctc)

# Tacks on links between the top nodes.
# Args:
#   clusters - the cluster labels (used to
#     determine names of clusters to add to)
#   r - the amount of correlation
#   f - the file to add on to
# Side effects: tacks a node onto the file
add.top.links = function(clusters, r, f) {
  num.clusters = max(clusters)

  # computes name of the i'th node to join
  ith.node = function(i)
    paste("NODE_", i, "_", sum(clusters==i)-1, "X", sep="")

  cat(paste("CLUSTER_NODE_1\t",
      ith.node(1), "\t", ith.node(2),
      "\t0\t#FFFFFF\n", sep=""), file=f, append=TRUE)

  if (num.clusters > 1) {
    for(i in 1:(max(clusters)-1)) {
      cat(paste("CLUSTER_NODE_", i+1, "\t",
        "CLUSTER_NODE_", i, "\t", ith.node(i+1),
          "\t0\t#FFFFFF\n", sep=""), file=f, append=TRUE)
    }
  }
}

# Writes out files.
# Args:
#   x - the expression dataset
#   wx - the expression dataset, weighted for clustering
#   clusters - which cluster each row is in
#   col - the color for each cluster
#   basefile - the base name of the output file
# Side effects: writes out the TreeView files
# Currently, this writes out a column clustering as well.
write.treeview = function(x, wx, clusters, col, basefile) {
  num.clusters = max(clusters)
  temp.file.base = "TreeViewTemp"

  for(s in c("atr", "gtr", "cdt")) {
    system(paste("rm ", basefile, ".", s, sep=""))
  }

  # header for gene tree
  cat("NODEID\tLEFT\tRIGHT\tCORRELATION\tNODECOLOR\n",
    file=paste(basefile, ".gtr", sep=""))

  # cluster arrays, and write this out
  # XXX for now, not sorting by this
  x1 = x
  x1[ is.na(x1) ] = 0
  hc <- hcluster(t(x1), method = "correlation", nbproc=4)
cat("computed hc\n")

  hc$order = sort(hc$order)
  r2atr(hc, file = paste(basefile, ".atr", sep = ""))

  # XXX without this, gene names are missing; I don't know why
  options(stringsAsFactors = TRUE)

  # cluster each cluster of genes
  for(i in 1:num.clusters) {

    wx1 = wx[ clusters == i , ]
    wx1[ is.na(wx1) ] = 0    # XXX hack

cat("i =", i, "\n")
cat("num genes =", sum(clusters==i), "\n")
cat("dim(wx1) =", dim(wx1), "\n")
    hr = hcluster(wx1, method = "correlation", nbproc = 7)

    # write out gene tree, colored appropriately
    r2gtr(hr, file = paste(temp.file.base, ".gtr", sep = ""))
    # tack that on to the output file with:
    # - NODE renamed per-cluster, and
    # - color tacked onto the end
    # XXX kind of a hack
    system(paste("perl -e '",
      "while (<>) {chomp; s/(NODE|GENE)/$1", "_", i, "_/g; ",
      "print \"$_\t", substr(col[i],1,7), "\n\"; }'", 
      "< ", temp.file.base, ".gtr >> ", basefile, ".gtr", sep=""))

#      "while (<>) {chomp; s/(NODE|GENE)/$1", "_", i, "_/; ",
#      "print \"$_\t\"", col[i], "\n\"; }'", 
#      "< ", temp.file.base, ".gtr > ", basefile, ".gtr", sep=""))

    # write out part of the data table
    r2cdt(hr, hc, x[ clusters==i,],
      file = paste(temp.file.base, ".cdt", sep = ""))

    # if this is the first cluster, copy the header
    if (i == 1) {
      system(paste("head -3 < ", temp.file.base, ".cdt",
        " >> ", basefile, ".cdt", sep=""))
    }

    # copy the other lines
    system(paste("perl -e '",
      "while (<>) {chomp; s/^GENE/GENE", "_", i, "_/; ",
      "if ($.>=4) { print \"$_\n\"; }; }'",
      "< ", temp.file.base, ".cdt >> ", basefile, ".cdt", sep=""))
  }

  add.top.links(clusters, -0.5, paste(basefile, ".gtr", sep=""))

  unlink(paste(temp.file.base, ".gtr", sep=""))
  unlink(paste(temp.file.base, ".ctd", sep=""))

  options(stringsAsFactors = FALSE)
}

# Alternative version of this, which just colors the original
# dendrogram. (This should put similar clusters nearby, making
# the overall clustering easier to interpret.) Somewhat hacky.
# Args:
#   x - the expression dataset to write out
#   hr, hc - row and column clustering to use
#   clusters - the cluster each gene is in
#   cluster.colors - the color for each cluster
#   basefile - base output file
#   label - the description to use for each gene
# Side effects: writes out files colored appropriately
write.clusters.treeview = function(x, hr, hc,
    clusters, cluster.colors, basefile, descr=NULL) {
  num.clusters = max(clusters)

  # first, write the files (not colorized)
  options(stringsAsFactors = TRUE)  # XXX don't know why this is needed
  r2gtr(hr, file = paste(basefile, ".gtr", sep = ""))
  r2atr(hc, file = paste(basefile, ".atr", sep = ""))
  if (is.null(descr))
    r2cdt(hr, hc, x, file = paste(basefile, ".cdt", sep = ""))
  else {
    x1 = cbind(rownames(x), as.character(descr),
      data.frame(x, check.names=FALSE))
# print("computed x1")
# print(class(x1))
# print(x1[1:5,1:5])
    r2cdt(hr, hc, x1, file = paste(basefile, ".cdt", sep = ""),
      labels=TRUE, description=TRUE)
  }

  # then, color the rows according to the clustering
  # read in files
  cdt = read.table(paste(basefile, ".cdt", sep = ""),
    sep="\t", quote="", skip=3, as.is=TRUE)
  gene.id = cdt[,1]
  names(gene.id) = cdt[,2]

  gtr = read.table(paste(basefile, ".gtr", sep = ""), as.is=TRUE)
  colnames(gtr) = c("NODEID", "LEFT", "RIGHT", "CORRELATION")
  rownames(gtr) = gtr$NODEID
  gtr$NODECOLOR = NA

  # vector for the names of nodes
  nn = unique(union(gtr$NODEID, union(gtr$LEFT, gtr$RIGHT)))
  node.colors = rep(NA, length(nn))
  names(node.colors) = nn

  # color genes according to their cluster
  # (note that TreeView can't deal with the alpha channel,
  # so these are truncated to, e.g., "#FF007A")
  node.colors[ gene.id[ names(clusters) ] ] =
    substr(cluster.colors[ clusters ], 1, 7)

  # propagate color up when a node's children are the same color
  done = FALSE
  while (!done) {
    nc1 = node.colors
    col.a = node.colors[ gtr[,"LEFT"] ]
    col.b = node.colors[ gtr[,"RIGHT"] ]
    i = col.a == col.b
    i[ is.na(i) ] = FALSE
    node.colors[ gtr[,"NODEID"][i] ] = col.a[ i ]   # FIXME
    gtr[,"NODECOLOR"] = node.colors[ gtr[,"NODEID"] ]
    done = all(node.colors == nc1, na.rm=TRUE) &&
      all(is.na(node.colors) == is.na(nc1))
  }
  gtr[,"NODECOLOR"] = node.colors[ gtr[,"NODEID"] ]

  # set remaining colors to grey
  gtr[ is.na(gtr[,"NODECOLOR"]), "NODECOLOR" ] = "#808080"

  # overwrite gene clustering with added color column
  write.table(gtr, file=paste(basefile, ".gtr", sep = ""),
    sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}

# Writes a "dummy" .atr file, which avoids changing the order
# of the samples.
# Args:
#   f - filename
#   n - number of "ARRY" entries (which will be numbered 0..n-1)
write.dummy.atr = function(f, n) {
  n = n - 1
  r = data.frame(NODEID=paste0("NODE", 1:n, "X"),
    LEFT = c("ARRY0X", paste0("NODE", 1:(n-1), "X")),
    RIGHT = paste0("ARRY", 1:n, "X"),
    CORRELATION = 1,
    NODECOLOR = "#ffffff")

  write.table(r, file=f, quote=FALSE, sep="\t",
    row.names=FALSE, col.names=TRUE)
}

