# Plots motif clusters.

library(grid)
library(MotIV)
library(latticeExtra)

source("git/plot/seqLogo.r")

load("git/tf/motif/clusterUsingMotIV.Rdata")

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
  lyt = grid.layout(n, 3,
    widths = unit(c(1,1,1), "in"),
    heights = unit(rep(0.3, n), "in"))
  pushViewport(viewport(layout=lyt))

  for(i in 1:n) {
    m = h$labels[ h$order[i] ] 

    # draw the labels
    pushViewport(viewport(layout.pos.row=i, layout.pos.col=1))
    grid.text(m, hjust=1)
    popViewport()

    # draw the motif logo
    pushViewport(viewport(unit(0, "in"), unit(0, "in"),
#      unit(1, "in"), unit(1, "in"),
      layout.pos.row=i, layout.pos.col=2))
    grid.rect(gp=gpar(col="lightgrey"))
    seqLogoSmall(pwm[[m]])
    popViewport()
  }

  # draw the dendrogram
  pushViewport(viewport(layout.pos.row=c(1:n), layout.pos.col=3))
  g = dendrogramGrob(as.dendrogram(h), side="right")
  grid.draw(g)
  popViewport()
  popViewport()
#  popViewport()
}

# Plots all the clusters.
plot.all = function() {
  pdf("git/tf/motif/plot/motif.clusters.pdf", width=7.5, height=10)

  m = c(jolma, de.novo.motifs)
  motif.groups = cutree(motif.clusters, k=50)
  for(i in 1:max(motif.groups)) {
cat(i, " ")

    motif.names = names(motif.groups)[motif.groups==i]
    m1 = list.select(m, motif.names)
    dists = motifDistances(m1)
    h = hclust(dists)

    grid.newpage()

    plot.motif.clustering(h, m)
  }

  dev.off()
}

plot.all()

