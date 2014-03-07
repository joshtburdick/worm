# Computes mean and s.d. of when genes are expressed.

source("git/utils.r")
source("git/data/name_convert.r")

# timeseries data (unscaled)
r.ts = read.tsv("git/unmix/seq/timing/deconvolved_embryo_ts.tsv")
r.ts = as.matrix(rename.gene.names(r.ts))

time.points = as.numeric(sub("t.", "", colnames(r.ts)))

# estimate mean and sd of when each gene is on
x1 = r.ts[ apply(r.ts, 1, mean) >= 1 , ]
w = x1 / apply(x1, 1, sum)
w[ is.na(w) ] = 0
gene.expr.time = {
  t.mean = apply( t(w) * time.points, 2, sum)
  t.s2 = apply( t(w) * (time.points^2), 2, sum)
  t.sd = sqrt( t.s2 - (t.mean^2) )
  data.frame(mean = t.mean, sd = t.sd)
}

gene.expr.time = gene.expr.time[ order(gene.expr.time$mean) , ]

write.tsv(gene.expr.time, "git/sort_paper/FACS/timing/geneExprTime.tsv")

