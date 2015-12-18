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

# reads per million (for computing average of maximum)
rpm = read.tsv("git/cluster/readsPerMillion.tsv")
rpm = rpm[ names(cl), !grepl("^(HS|RNAi)", colnames(rpm)) ]

# the Spencer data
load("data/expression/spencer.expr.Rdata")
s1 = spencer.expr[[1]]
rownames(s1) = s1$name
s1 = as.matrix(s1[,4:33])
colnames(s1)[4] = "LE.reference"
s1 = rename.gene.names(s1)

# limit to genes in common
g = intersect(rownames(s1), names(cl))
s1 = s1[ g , ]
r.facs.1 = r.facs[ g , ]
cl.s = cl[ g ]
# the Spencer data, normalized
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

# average max. expression of each gene
max.rpm = apply(rpm, 1, max)
mean.max.lrpm = tapply(log2(1+max.rpm), cl, mean)
mean.max.lrpm = mean.max.lrpm[ order(as.numeric(names(mean.max.lrpm))) ]
max.spencer = apply(s1, 1, max)
mean.max.spencer = tapply(max.spencer, se.cl, mean)

cluster.coloring = rep(hsv(c(1:5) / 6, 1, 0.7), 100)

# Compares correlations between clusterings.
compare.correlation = function(a.cl, a, b, basefile, a.name, b.name, cl.expression, cl.to.color) {
  a = a[ names(a.cl) , ]
  b = b[ names(a.cl) , ]
# browser()
  mc.clustered = by(a, a.cl, mean.cor)
  mc.other = by(b, a.cl, mean.cor)
  rlim = c(-0.1, 1)

  # coloring for clusters
  expr.range = as.numeric(round(quantile(cl.expression, c(0.05, 0.95))))
  cl.expression.1 = (cl.expression - expr.range[1]) / (expr.range[2] - expr.range[1])
  cl.expression.1[ cl.expression.1 < 0 ] = 0
  cl.expression.1[ cl.expression.1 > 1 ] = 1
  cl.expression.1 = cl.expression.1[ names(mc.clustered) ]

  col.f = function(x) hsv(0, 0, 0.8 * x, 0.8)
  cl.col = col.f(1 - cl.expression.1)
  names(cl.col) = names(cl.expression.1)
  cl.col[ cl.to.color ] = hsv(0, 1, 1, 0.8)

  plot(mc.clustered, mc.other, xlim=rlim, ylim=rlim, pch=20,
    col = cl.col, cex=0.8,
    xlab=paste0("Correlation in ", a.name, " data"),
    ylab=paste0("Correlation in ", b.name, " data"))

  # legend for expression scale
  x.legend = c(expr.range[2]:expr.range[1])
  x.legend.scaled = (x.legend - expr.range[1]) / (expr.range[2] - expr.range[1])
  legend(-0.1, 0.98, legend = x.legend, pch = 20, col = col.f(x.legend.scaled),
    cex=0.7, bty="n")
  text(-0.1,1, adj=c(0,0.5), cex=0.7, ifelse(a.name=="Spencer",
      expression("Average expression (max. intensity)"),
      expression("Average expression (log"[2] * "(1 + max. RPM))")))

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



if (TRUE) {
pdf("git/sort_paper/cluster/spencerComparison/clusterComparison.pdf",
  width=8, height=4)
par(mar=c(5,5,1,1))
par(mfrow=c(1,2))
by.spencer = compare.correlation(se.cl, se, r.facs.1, paste0(base.dir, "bySpencer"),
  "Spencer", "FACS", mean.max.spencer, "244")
label.panel("A")
by.facs = compare.correlation(cl.s, r.facs.1, se, paste0(base.dir, "byFacs"),
  "FACS", "Spencer", mean.max.lrpm, "189")
label.panel("B")

dev.off()
}

