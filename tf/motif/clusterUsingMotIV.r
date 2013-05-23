# Clusters motifs using MotIV.

library("MotIV")
library("seqLogo")

source("git/plot/seqLogo.r")

# Converts a list of motifs to a matrix, padded with NAs.
# Args:
#   m - a list of motifs (as returned by, e.g., readPWMfile),
#     each entry of which is a 4-by-w matrix (where w is
#     the motif length), with rownames A, C, G, and T.
# Returns: a 4-by-maxWidth-by-n array, in which
#   maxWidth is the width of the widest motif, and
#   n is the number of sequences.
# (Not yet implemented.)



# hack to read a MEME format file directly.
read.meme.file = function(f) {
  # XXX should use a proper temp file, and, like, check for errors

#  system(paste("/home/jburdick/meme/bin/meme2meme ", f, " > /var/tmp/meme.txt"))
#  system(paste("/home/jburdick/gcb/git/perl/motif/meme2TRANSFAC.pl /var/tmp/meme.txt /var/tmp/meme.transfac.txt")
  system(paste("/home/jburdick/gcb/git/perl/motif/meme2TRANSFAC.pl ", f, " /var/tmp/meme.transfac.txt"))
  r = readPWMfile("/var/tmp/meme.transfac.txt")
  r
}

# ??? write a function to read MEME-format files directly?
# a = readPWMfile("git/perl/motif/hier.ts.200clusters.TRANSFAC.txt")

a = read.meme.file("git/tf/motif/meme_file/hier.ts.200clusters.meme")
jolma = read.meme.file("data/tf/meme/motif_databases/jolma2013.meme")

# Selects some things from a list.
# ??? is there a function which would do this?
list.select = function(a, x) {
  r = lapply(x, function(x) a[[x]])
  names(r) = x
  r
}

# Reduces a set of clusters using a particular cutoff.
cluster.similar.motifs = function(motifs, height) {
  dists = motifDistances(motifs)
cat("computed distances\n")
  h = hclust(dists)
  motif.groups = cutree(h, h=height)
  list.select(motifs, names(motif.groups)[!duplicated(motif.groups)])
}

# Plots groups of motifs which are within some similarity cutoff.
# XXX somewhat deprecated
plot.motif.clusters = function(m, num.clusters, output.file) {

  dists = motifDistances(m)
  h = hclust(dists)
  motif.groups = cutree(h, k = num.clusters)

  pdf(output.file, width=11, height=8.5)


  for (i in 1:5) {    # max(motif.groups)) {
#    grid.newpage()
  pushViewport(viewport(0.2, 0.2, 0.6, 0.6,
    default.units="npc", clip="off"))
cat(i, " ")
    motif.names = names(motif.groups)[motif.groups==i]
    if (length(motif.names) > 1) {
      m1 = list.select(m, motif.names)
      dists = motifDistances(m1)
      h = hclust(dists)
# rotate this
#vp = viewport(0.5, 0.5, 1, 1, angle=90, default.units="npc",
#  just=c("center", "center"), clip="off")
#pushViewport(vp)
  dendro = as.dendrogram(h)
    par(mar=c(6,4,4,8)+0.2)
      plot(dendro, cex=0.7, horiz=TRUE)
      draw.hclust.seq.logos(h, m1)
#      popViewport()
    }
    else {
      cat("skipping ")
    }
    popViewport()
  }
  dev.off()
}

# quick test
# seqLogoSmall(a[[1]])


a1 = cluster.similar.motifs(a, 0.1)
m = c(jolma, a1)
plot.motif.clusters(m, 40, "git/tf/motif/jolmaClustered.pdf")





if (FALSE) {


print(system.time(dists <- motifDistances(a)))

# cluster motifs, and pick one from each cluster
h = hclust(dists)
a1 = list()
motif.groups = cutree(h, h=0.001)
for(n in names(motif.groups)[!duplicated(motif.groups)]) {
  a1[[n]] = a[[n]]
}

}

# FIXME use these?
draw.clusters.logo = function() {
  grid.newpage()
  vp <- viewport(x=.5, y=.5,
    width=.5, height=.5,
    just=c("left", "bottom"))
  pushViewport(vp)
#  grid.rect(gp=gpar(fill="grey"))

#  plot(a1)
   popViewport(1)
}

