# Unmixing including temporal data.

source("git/unmix/unmix.r")
source("git/plot/plot_expr.r")

embryo.timeseries = read.table("git/unmix/seq/timing/embryo.timeseries.tsv.gz",
  sep="\t", header=TRUE, row.names=1, as.is=TRUE)
time.points = c(0,30,60,90,120,140,180,240,270,300,330,360,390,
  420,450,480,540,570,600,630,660,690,720)
# colnames(embryo.timeseries) = time.points
g1 = intersect(rownames(embryo.timeseries), rownames(r.corrected$r.mean))
embryo.timeseries = embryo.timeseries[ g1 , ]

# normalize to "ppm"
embryo.timeseries =
  t( t(embryo.timeseries) / (apply(embryo.timeseries,2,sum) / 1e6) )

x.ft = cbind(r.corrected$r.mean[g1,], embryo.timeseries[g1, c(2:18)])

load("git/unmix/image/time_sort_matrix.Rdata")

# scale rows of this to add up to 1
M.t = time.sort.matrix / apply(time.sort.matrix, 1, sum)
M.t = rbind(M, M.t)

x.p.time = x.ft %*% pseudoinverse(t(M.t))

x.p.time[ is.na(x.p.time) ] = 0
x.p.time[ x.p.time < 0 ] = 0

plot.it = function(genes) {
  pdf("git/unmix/unmix_time.pdf", width=10, height=10)
  par(mfrow=c(2,1))

  for(g in genes) {
    cat(g, "")
    r = rgb(scale.to.unit(x.pseudoinverse[g,]), 0, 0)
    names(r) = colnames(x.pseudoinverse)
    plot.segments.per.cell(r, paste(g, "(pseudoinverse)"),
      times=c(-30, 580), lwd=3)
  
    r = rgb(scale.to.unit(x.p.time[g,]), 0, 0)
    names(r) = colnames(x.p.time)
    plot.segments.per.cell(r, paste(g, "(pseudoinverse with time)"),
      times=c(-30, 580), lwd=3)
  }

  dev.off()
}

# get list of genes to include
load("R/unmix/comp_paper/expr.cell.Rdata")

enriched.fraction = read.table("R/unmix/sort_paper/unmix/fraction/enriched.fraction.tsv", sep="\t", header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE)

gene.list = c(rownames(expr.cell), rownames(enriched.fraction))
gene.list = unique(intersect(gene.list, rownames(x.ft)))

x.pc.time = unmix.lsei(M.t, x.ft[gene.list,], x.ft[gene.list,])


# plot.it(sample(rownames(enriched.fraction), 5))

