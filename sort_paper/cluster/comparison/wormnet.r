# Comparison of clustering with Wormnet annotation.

source("git/utils.r")
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

# the read ratios
x = read.tsv("git/cluster/readRatios.tsv")[,1:23]

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

# this was:
#  p.cl = table(cl)
#  p.cl = p.cl / sum(p.cl)
#  background.prob = sum(p.cl^2)

  # but I think it should be the following (as it only counts
  # a given interaction once)
  n.cl = table(cl)
  background.prob = sum(choose(n.cl, 2)) / (sum(n.cl)^2)

  for(arc.type in arc.types) {
    s = w1[ !is.na(w1[,arc.type]) , "same.cluster" ]
    a = chisq.test(c(sum(s), sum(!s)),
      p=c(background.prob, 1-background.prob))
    r = rbind(r,
      data.frame(
        num.arcs = sum(!is.na(w1[,arc.type])),
        num.same.cluster = sum( s ),
        p.same.cluster =
          signif(mean(w1[ !is.na(w1[,arc.type]), "same.cluster" ]), 2),
        chisq = as.numeric(a$statistic),
        p = a$p.value
     ))
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

# Alternative version of this, which compares the number
# overlapping with a background based on shuffling
# (a la Walhout et al 2013).
# Args:
#   a - arcs, as a data frame with two columns
#   cl - a clustering
#   num.shuffles - number of times to sample this
# Returns: vector with elements
#   n.same.cluster - number of arcs in the same cluster 
#   shuffled.mean - mean of that, for all shuffles
#   n.shuffled.higher
cluster.overlap.shuffle.1 = function(a, cl, num.shuffles=100) {
  n.same.cluster = sum(cl[a[,1]] == cl[a[,2]])
  r = rep(NA, num.shuffles)

  for(i in 1:num.shuffles) {
    if (i %% 100 == 0) {
      write.status(i)
    }
    s = sample(cl)
    names(s) = names(cl)
    r[i] = sum(s[a[,1]] == s[a[,2]])
  }
  n.shuffled.higher = sum(r >= n.same.cluster)

  data.frame(n.same.cluster = n.same.cluster,
#    n.same.cluster.shuffled = r,
    shuffled.mean = mean(r),
    n.shuffled.higher = n.shuffled.higher,
    p = (n.shuffled.higher + 1) / num.shuffles,
    stringsAsFactors=FALSE)
}

# Checks the cluster overlap for each type of arc.
cluster.overlap.shuffle = function(num.shuffles=10) {
  r = NULL

  for(arc.type in arc.types) {
    cat(arc.type, "\n")
    s = w1[ !is.na(w1[,arc.type]) , c("gene1","gene2") ]
    r1 = cluster.overlap.shuffle.1(s, cl, num.shuffles)
    r = rbind(r, r1)
  }
  rownames(r) = arc.types
  r$p.corr = p.adjust(r$p, method="fdr")
  r
}

# original version of this
if (FALSE) {
arc.cluster.counts = compare.arc.clusters()

arc.cluster.counts[ 1:21 , ] =
  arc.cluster.counts[
    order(arc.cluster.counts[1:21,2], decreasing=TRUE) , ]

write.tsv(arc.cluster.counts, "git/sort_paper/cluster/comparison/wormnet.tsv")
}

if (FALSE) {
 z = w1[ !is.na(w1$"SC-CX") , ]
z1 = cluster.overlap.shuffle.1(z, cl, 20000)

}

wormnet.stats.shuffle = cluster.overlap.shuffle(10000)

write.tsv(wormnet.stats.shuffle,
  "git/sort_paper/cluster/comparison/wormnet_shuffle.tsv")

