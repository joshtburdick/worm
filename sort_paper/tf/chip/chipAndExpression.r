# Compares average expression and enrichment of
# ChIP motifs.

source("git/utils.r")

# directory containing the ChIP data
chip.count.dir = "git/tf/chip/count/upstreamChipCount/"

# upstream region sizes
load("git/tf/motif/count/upstreamRegionSize.Rdata")

# the read ratios
x = as.matrix(read.tsv("git/cluster/readsPerMillion.tsv"))
x = x[ , grep("HS |RNAi", colnames(x), invert=TRUE) ]

# Computes center of each cluster.
# Args:
#   r - the expression data
#   cl - the clustering (as a vector of integers,
#     with genes as names)
# Returns: center of each cluster
compute.cluster.centroids = function(r, cl) {
  cluster.names = sort(unique(cl))
  num.clusters = max(cl)

  # first, compute cluster centers
  r.center = matrix(nrow=num.clusters, ncol=ncol(r))
  rownames(r.center) = cluster.names
  colnames(r.center) = colnames(r)

  for (i in cluster.names)
    r.center[i,] = apply(r[cl==i,], 2,
      function(x) mean(x, na.rm=TRUE))

  r.center
}

clustering = "hier.300.clusters"
cl1 = read.tsv(paste0("git/cluster/hierarchical/", clustering, "/clusters.tsv"))
cl = cl1[,2]
names(cl) = rownames(cl1)
x1 = x[names(cl),]

# compute log-transformed mean expression
lme = log2(1 + compute.cluster.centroids(x1, cl))

cluster.lme = apply(lme, 1, mean)

# Totals motifs and upstream region sizes in clusters.
count.in.clusters = function(cl, motif.count, upstream.bp) {
  motif.count = motif.count[names(cl),,,,drop=FALSE]
  upstream.bp = upstream.bp[names(cl),,]

# browser()

  # count motifs by cluster (foreground and background)
  motifs.cluster =
    aperm(apply(motif.count, c(2,3,4), function(x) by(x, cl, sum)), c(2,3,4,1))
  motifs.total = apply(motifs.cluster, c(1,2,3), sum)
  # aperm() was used so that this works
#  motifs.background = as.vector(motifs.total) - motifs.cluster

  # add up amount of upstream sequence (foreground and background)
  bp.cluster = aperm(apply(upstream.bp, c(2,3),
    function(x) by(x, cl, sum)), c(2,3,1))
  bp.total = apply(bp.cluster, c(1,2), sum)
#  bp.background = as.vector(bp.total) - bp.cluster

  # more aperm() gymnastics
  motifs.cluster = aperm(motifs.cluster, c(4,1,2,3))
#  motifs.background = aperm(motifs.background, c(4,1,2,3))
  bp.cluster = aperm(bp.cluster, c(3,1,2))
#  bp.background = aperm(bp.background, c(3,1,2))

  list(motifs.cluster = motifs.cluster, bp.cluster = bp.cluster)
}

# Correlates one ChIP experiment with expression.
chip.expr.cor = function(chip.count, bp.count, cluster.lme) {
  r = NULL

  # XXX not sure why this is happening
  if (sum(chip.count)==0)
    return(NULL)

  for(ud in dimnames(chip.count)[[2]])
    for(cons in dimnames(chip.count)[[3]])
      for(score in dimnames(chip.count)[[4]]) {
        peaks.per.kb = chip.count[,ud,cons,1] /
          (bp.count[,ud,cons] / 1e3)
        if (any(peaks.per.kb > 0)) {
          s = summary(lm(cluster.lme ~ peaks.per.kb))
          r = rbind(r, data.frame(upstream.dist = ud,
            conservation = cons,
            t = s$coefficients[2,3],
            p = s$coefficients[2,4],
            r.squared = s$r.squared))
        }
      }
  r
}

# Does this correlation for all the ChIP experiments.
chip.expr.cor.all = function(chip.dir, cluster.lme) {
  r = NULL

  for(f in list.files(chip.dir)) {
    chip.experiment = sub(".Rdata$", "", f)
    write.status(chip.experiment)
    motif.count = NULL
    load(paste0(chip.count.dir, f))
    n = count.in.clusters(cl, motif.count, upstream.region.size)
    r1 = cbind("chip.experiment"=chip.experiment,
      chip.expr.cor(n$motifs.cluster, n$bp.cluster, cluster.lme),
      stringsAsFactors=FALSE)
    r = rbind(r, r1)
  }
  cat("\n")

  r
}

r = chip.expr.cor.all(chip.count.dir, cluster.lme)

# adjust p-values
r$p.corr = p.adjust(r$p, method="fdr")

# just keep the most significant for each experiment
r = r[ order(r$p.corr) , ]
r = r[ !duplicated(r[,c("chip.experiment")]) , ]

write.tsv(r, "git/sort_paper/tf/chip/chipAndExpression.tsv")

