# Merges WGCNA dendrograms.
# XXX not working; possibly won't be used.

library(WGCNA)

# Merges two hclust objects.
# Args: a, b - two hclust objects; these should have
#   non-overlapping labels.
# Returns: an hclust object, with these merged
# (the height of the joining node is 1)
merge.hclust = function(a, b) {
  # number of leaves in each
  m.a = length(a$order)
  m.b = length(b$order)

  # shift numbers in b$merge
  m1 = b$merge
  m1[ m1 < 0 ] = m1[ m1 < 0 ] - m.a
  m1[ m1 > 0 ] = m1[ m1 > 0 ] + (m.a - 1)

  r = list(
    merge = rbind(a$merge, m1, c(m.a - 1, m.a + m.b - 2)),
    height = c(a$height, b$height, 1),
    order = c(a$order, b$order + m.a),
    labels = c(a$labels, b$labels),
    method = a$method, call = a$call, distance = a$distance)
  class(r) = "hclust"
  r
}

# Sorts an hclust object so that heights are in ascending order.
sort.hclust = function(h) {
  p = order(h$height)

  # renumber the nodes
  m1 = h$merge
  m1[ m1 > 0 ] = order(p)[ m1[ m1 > 0 ] ]

  r = list(
    merge = m1[ p , ],
    height = h$height[ p ],
    order = h$order,
    labels = h$labels,
    method = h$method, call = h$call, distance = h$distance)
  class(r) = "hclust"
  r
}

# Generates a random hclust (for testing.)
random.hclust = function(n) {
  x = runif(n)
  names(x) = paste("x", sample(1000, n), sep="")
  hclust(dist(x))
}

# Quick test of the above.
merge.hclust.test = function(n) {

  h1 = random.hclust(n)
  h2 = random.hclust(n)
  h3 = merge.hclust(h1, h2)
  h3s = sort.hclust(h3)
  par(mfrow=c(2,2))
  plot(h1, ylim=c(0,1))
  plot(h2, ylim=c(0,1))
  plot(h3, ylim=c(0,1))
  plot(h3s, ylim=c(0,1))
  list(h1=h1, h2=h2, h3=h3, h3s=h3s)  
}


# Merges dendrograms (deprecated.)
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


