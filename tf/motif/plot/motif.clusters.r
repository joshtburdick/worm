# Plots motif clusters.
# This is fairly inefficient code, as it does a lot of
# re-clustering; hopefully the code is thus slightly more
# straightforward.

library(grid)
library(MotIV)
library(latticeExtra)

source("git/plot/seqLogo.r")

# read in the clustering
# load("git/tf/motif/clusterUsingMotIV.Rdata")
# load("git/tf/motif/motif.comparison.results.Rdata")
load("git/tf/motif/clusterUsingMotifComparison.Rdata")

# read in the motifs
load("git/tf/motif/meme.format.pwm.Rdata")

# Gets a subset of the rows from the MotifComparison results,
# and clusters them hierarchically.
# Args:
#   m - a set of motif names
# Returns: an hclust object from grouping them
motif.subset.hclust = function(m) {
  r = motif.comparison.results
  r1 = r[ r$a %in% m & r$b %in% m , ]
  hclust(MCFile.to.dist(r1))
}

# Converts a MotifComparison output file to a distance matrix.
MCFile.to.dist = function(r) {
  motif.names = sort(unique(c(r$a, r$b)))
  a = matrix(4, nrow=length(motif.names), ncol=length(motif.names))
  rownames(a) = motif.names
  colnames(a) = motif.names

  a[cbind(r$a, r$b)] = r$score
  a[cbind(r$b, r$a)] = r$score
  as.dist(a)
}

# Utility to select things from a list.
# ??? is there a function which would do this?
list.select = function(a, x) {
  r = lapply(x, function(x) a[[x]])
  names(r) = x
  r
}

# Plots a clustering of motifs.
# Args:
#   h - the hierarchical clustering to draw
#   pwm - a list including PWMs for these motifs (and possibly
#     others), whose names should agree with h$labels
# Side effects: draws a diagram with names of motifs,
#   motif logos, and a dendrogram
plot.motif.clustering = function(h, pwm) {
  n = length(h$order)

#  pushViewport(viewport(unit(5, "in"), unit(0.2, "in")))
  lyt = grid.layout(n+2, 3,
    widths = unit(c(1.5,1,1.5), "in"),
    heights = unit(rep(0.2, n+2), "in"))
  pushViewport(viewport(layout=lyt))

  for(i in 1:n) {
    m = h$labels[ h$order[i] ]     # XXX
    # draw the labels
    pushViewport(viewport(layout.pos.row=(n+1)-i, layout.pos.col=1))
    grid.text(m, hjust=1, gp=gpar(cex=0.95))
    popViewport()

    # draw the motif logo
    pushViewport(viewport(unit(0, "in"), unit(0, "in"),
#      unit(1, "in"), unit(1, "in"),
      layout.pos.row=(n+1)-i, layout.pos.col=2))
    grid.rect(gp=gpar(col="lightgrey"))
    seqLogoSmall(pwm[[m]])
    popViewport()
  }

  # draw the dendrogram
  pushViewport(viewport(layout.pos.row=c(1:n), layout.pos.col=3))
  dendro = as.dendrogram(h)
  dg = dendrogramGrob(dendro, side="right")
  grid.draw(dg)
#  grid.xaxis(main=TRUE, gp=gpar(cex=0.7))
  popViewport()

  # draw the scale for the dendrogram
  # FIXME: this isn't aligned correctly...
if (FALSE) {
  pushViewport(viewport(layout.pos.row=c(n), layout.pos.col=3,
    xscale=c(0,max(h$height))))
  grid.xaxis(main=TRUE, gp=gpar(cex=0.7))
  popViewport()
  pushViewport(viewport(layout.pos.row=c(n+2), layout.pos.col=3))
  grid.text("average KL distance",
    x=unit(0.5, "npc"), y=unit(0, "npc"), gp=gpar(cex=0.8))
  popViewport()
}
# for testing. XXX throwing an exception, dunno why...
  if (n > 2) {
    try({
      plot(h)
    })
  }

#  popViewport()
}

# Plots all the clusters.
plot.all = function() {
  pdf("git/tf/motif/plot/motif.clusters.pdf", width=7.5, height=10)

#  m = cluster.using.motif.comparison$labels
  motif.groups = cutree(cluster.using.motif.comparison, k=300)
  for(i in 1:max(motif.groups)) {
cat(i, " ")

    motif.names = names(motif.groups)[motif.groups==i]
    h = motif.subset.hclust(motif.names)

    grid.newpage()

    plot.motif.clustering(h, meme.format.pwm)
  }

  dev.off()
}

plot.all()

