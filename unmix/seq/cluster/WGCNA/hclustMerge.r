# Merges WGCNA dendrograms.
# XXX not working; possibly won't be used.

library(WGCNA)

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

# Merges dendrograms from a WGCNA analysis.
# Args:
#   wnet - the list returned by blockwiseModules
#   name - the names of the genes
# Returns: a dendrogram
get.wgcna.dendrogram = function(wnet, name) {
  numBlocks = length(wnet$blockGenes)

  # add labels to hclust objects
  h = list()
  for(i in 1:numBlocks) {
    h[[i]] = wnet$dendrograms[[i]]
    h[[i]]$labels = name[ wnet$blockGenes[[i]] ]
  }

  # if there's only one block, we're done
  if (numBlocks == 1)
    return(h[[1]])

  # otherwise, merge the hclust objects together
  r = h[[1]]
  for(i in 2:numBlocks)
    r = merge.hclust(r, h[[i]])

  sort.hclust(r)
}

