# Another comparison of clustering.
# in Spencer, cluster 244 with flp-13
# in FACS, cluster 84 with col-98?
#   cluster 189 with asm-1?

source("git/utils.r")
source("git/plot/label_panel.r")
source("git/data/name_convert.r")

source("git/unmix/seq/cluster/writeClustersTreeView.r")

output.path = "git/sort_paper/cluster/spencerComparison/clusterComparison"

# read ratios, from FACS data
r = as.matrix(read.tsv("git/cluster/readRatios.tsv"))

# clustering
clustering = "hier.300.clusters"
cl1 = read.tsv(paste0("git/cluster/hierarchical/", clustering, "/clusters.tsv"))
cl = cl1[,2]
names(cl) = rownames(cl1)
r1 = r[names(cl),]

r.facs = r1[ , c(1:23) ]

# the Spencer data
load("data/expression/spencer.expr.Rdata")
s1 = spencer.expr[[1]]
rownames(s1) = s1$name
s1 = as.matrix(s1[,4:33])
colnames(s1)[4] = "LE.reference"
s1 = rename.gene.names(s1)

g = intersect(rownames(s1), names(cl))
s1 = s1[ g , ]
r.facs.1 = r.facs[ g , ]
cl.s = cl[ g ]

se = cbind( s1[,c(2:3)] - s1[,1], s1[,c(5:15)] - s1[,4] )

# do clustering
se.hclust = hcluster(se, method="correlation", link="complete", nbproc=4)
se.cl = cutree(se.hclust, k=300)

# the mean correlation between distinct rows of a matrix
mean.cor = function(x) {
  r = cor(t(x))
  diag(r) = NA
  mean(as.vector(r), na.rm=TRUE)
}

cluster.coloring = rep(hsv(c(1:5) / 6, 1, 0.7), 100)

# Compares correlations between clusterings.
compare.correlation = function(a.cl, a, b, basefile, a.name, b.name) {
  a = a[ names(a.cl) , ]
  b = b[ names(a.cl) , ]
# browser()
  mc.clustered = by(a, a.cl, mean.cor)
  mc.other = by(b, a.cl, mean.cor)
  rlim = c(-0.1, 1)
  plot(mc.clustered, mc.other, xlim=rlim, ylim=rlim, pch=20, col="#50505080",
    xlab=paste0("Correlation in ", a.name, " data"),
    ylab=paste0("Correlation in ", b.name, " data"))

  # write out Treeview files (in a given order)
  combined.data = cbind(r.facs.1, " "=NA, se)
#  combined.data = cbind(a, " "=NA, b)
  write.treeview(combined.data, a, a.cl, cluster.coloring, basefile)

  # add on a representative gene per cluster (adding cluster labels to
  # the file might be better)
  sample.gene = c(by(names(a.cl), a.cl, function(x) x[1]))
  r = cbind(clustered = mc.clustered, other=mc.other)
  r = as.data.frame(r)
  r$sample.gene = sample.gene[rownames(r)]
  r = r[ order( r$clustered - r$other, decreasing = TRUE ) , ]
  r
}

base.dir = "git/sort_paper/cluster/spencerComparison/"
system(paste0("mkdir -p ", base.dir))


pdf("git/sort_paper/cluster/spencerComparison/clusterComparison.pdf",
  width=8, height=4)
par(mar=c(5,5,1,1))
par(mfrow=c(1,2))
by.spencer = compare.correlation(se.cl, se, r.facs.1, paste0(base.dir, "bySpencer"),
  "Spencer", "FACS")
label.panel("A")
by.facs = compare.correlation(cl.s, r.facs.1, se, paste0(base.dir, "byFacs"),
  "FACS", "Spencer")
label.panel("B")

dev.off()


