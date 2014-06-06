# Comparison of clustering with Wormnet annotation.

source("git/data/worm/wormnet.r")

arc.types = colnames(wormnet)[3:24]

# one of the clusterings
cl = {
  cl1 = read.table("git/cluster/hierarchical/hier.300.clusters/clusters.tsv",
    sep="\t", header=TRUE, as.is=TRUE)
  cl = cl1[,3]
  names(cl) = cl1[,2]
  cl
}

# only keep the overlapping gene names (the overlap seems small,
# but then, the clustering doesn't include genes not expressed
# in the embryo, for instance)
g = intersect(names(cl), union(wormnet$gene1, wormnet$gene2))
w1 = wormnet[ (wormnet$gene1 %in% g) & (wormnet$gene2 %in% g) , ]
w1$same.cluster = (cl[w1$gene1] == cl[w1$gene2])
cl = cl[g]

# Plots a comparison of the scores for these; they don't seem
# very distinct.
plot.score.dist = function() {
  pdf("git/sort_paper/cluster/comparison/wormnet.pdf",
    width=10, height=7.5)
  par(mfcol=c(2,4))

  for(arc.type in colnames(w1)[3:24]) {
    lim = range(w1[ , arc.type ])
    hist(w1[ ! w1$same.cluster, arc.type ], col="lightgrey", breaks=50,
      main = "Different clusters", xlab=arc.type)
    hist(w1[ w1$same.cluster, arc.type ], col="darkgrey", breaks=50,
      main = "Same cluster", xlab=arc.type)
  }

  dev.off()
}

# Compares how often ends of an arc are in the same cluster.
compare.arc.clusters = function() {
  r = NULL

  p.cl = table(cl)
  p.cl = p.cl / sum(p.cl)
  background.prob = sum(p.cl^2)

  for(arc.type in arc.types) {
    r = rbind(r,
      data.frame(
        p.same.cluster = round(mean(w1[ !is.na(w1[,arc.type]), "same.cluster" ]), 4) ))
  }

  rownames(r) = arc.types
  r$times.more.than.random = round(r$p.same.cluster / background.prob, 4)
  r$description = wormnet.evidence.type[ arc.types ]
  r
}






# Comparison of enrichment of genetic interactions.
compare.genetic.interactions = function() {

  genes = unique(c(w1.b[,1], w1.b[,2]))
  
  p.cx = sum(w1.b$"CE-CX") / choose(length(genes), 2)
  cl1 = cl[genes]
  p.cl = table(cl1)
  p.cl = p.cl / sum(p.cl)
  p.same.cluster = sum(p.cl^2)
  n2 = choose(length(genes), 2)
  num.same.cluster = sum(choose(table(cl1), 2))

  list(num.genes = length(genes),
    p.ce.gn = sum(w1.b[,"CE-GN"]) / choose(length(genes),2),
    p.ce.gn.given.cx = sum(w1.b[ w1.b[,"CE-CX"],"CE-GN"]) / sum(w1.b[,"CE-CX"]),
    p.ce.gn.given.same.cluster = sum(w1.b[ w1.b[,"same.cluster"],"CE-GN"]) / num.same.cluster)
}

arc.cluster.counts = compare.arc.clusters()
write.tsv(arc.cluster.counts, "git/sort_paper/cluster/comparison/wormnet.tsv")

w1.b = w1
w1.b[ , c(3:24) ] = !is.na( w1.b[ , c(3:24) ] )





