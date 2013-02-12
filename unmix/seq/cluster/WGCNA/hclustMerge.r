# Merges WGCNA dendrograms.
# XXX not working; possibly won't be used.

library(WGCNA)

# Merges dendrograms.
mergeWGCNAdendrograms = function(wnet) {
  numBlocks = length(wnet$blockGenes)

  # number of genes in each block
  blockSizes = sapply(wnet$blockGenes, length)

  # amount to add to, or actually subtract from, leaves
  leafOffsets = c(0, cumsum(blockSizes[-numBlocks]))
  
  # amount to add to internal nodes
  internalOffsets = c(0, cumsum(blockSizes[-numBlocks]-1))

  r = list(merge = c(),
    height = NULL,
    order = NULL,
    labels = NULL,
    call = wnet$dendrograms[[1]]$call,
    method = wnet$dendrograms[[1]]$method)

  for(i in 1) {    # c(1:numBlocks)) {

    m1 = wnet$dendrograms[[i]]$merge

    # renumber leaf nodes
    m1[ m1 < 0 ] = m1[ m1 < 0 ] - leafOffsets[i]
#    m1[m1 < 0] = - wnet$blockGenes[[i]][
#      wnet$dendrograms[[i]]$order[ - m1[m1 < 0] ] ]

    # renumber internal nodes
    m1[m1 > 0] = m1[m1 > 0] + internalOffsets[i]

    r$merge = rbind(r$merge, m1)
    r$height = c(r$height, wnet$dendrograms[[i]]$height)
    r$order = c(r$order, wnet$blockGenes[[i]])
#    r$order = c(r$order, wnet$blockGenes[[i]][ wnet$dendrograms[[i]]$order ])  #+ blockOffsets[i])
  }

  # add in any genes that haven't been added
  r$order = c(r$order, setdiff(1:100, r$order))
#  r$order = c(1:100)

  # FIXME add links between all clusters

  class(r) = "hclust"
  r
}



