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
  hc <- hcluster(t(x), method = "correlation", nbproc=4)
  r2atr(hc, file = paste(basefile, ".atr", sep = ""))

  options(stringsAsFactors = TRUE)
  # cluster each cluster of genes
  for(i in 1:num.clusters) {

    wx1 = wx[ clusters == i , ]

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


