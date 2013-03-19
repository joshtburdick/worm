# Using the deconvolution of the embryonic timeseries.

library("ctc")
library("biomaRt")

wb.gene = useMart("WS220", "wormbase_gene")

# Generalized "conversion" function.
# Args:
#   mart - a biomaRt object
#   from - attribute to filter on
#   to - attribute to return
# Returns: function which converts names.
mart.convert = function(mart, from, to) function(x) {
  r = getBM(attributes = c(from, to),
    filters = c(from),   
    values = list(x),
    mart = mart)
  r
  r1 = r[ match(x, r[,1]) , 2 ]
  r1
}

# original quantification, from .wig files
embryo.timeseries = as.matrix(read.table(
  "git/unmix/seq/timing/embryo.timeseries.tsv.gz",
  sep="\t", header=TRUE, row.names=1, as.is=TRUE))
embryo.timeseries =
  t( t(embryo.timeseries) / (apply(embryo.timeseries, 2, sum) / 1e6) )


# read counts from Tophat mapping
r.tophat = as.matrix(
  read.table("git/unmix/seq/quant/readsPerMillion/readsPerMillion_embryo_ts.tsv"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE)
colnames(r.tophat) = sub("50.0_", "50.000_", colnames(r.tophat))
colnames(r.tophat) = sub("50.30_", "50.030_", colnames(r.tophat))
colnames(r.tophat) = sub("50.60_", "50.060_", colnames(r.tophat))
colnames(r.tophat) = sub("50.90_", "50.090_", colnames(r.tophat))
r.tophat = r.tophat[ , sort(colnames(r.tophat)) ]

# deconvolution from modENCODE
r.deconvolve = as.matrix(read.table(gzfile("data/modencode/N2_EE_50_mod1.csv.gz"),
  sep=",", header=FALSE, row.names=1))
colnames(r.deconvolve) = c("t.000", "t.060", "t.120", "t.150", "t.180",
  "t.240", "t.330", "t.390", "t.420", "t.480",
  "t.540", "t.570", "t.600", "t.630", "t.660")

rownames(r.deconvolve) =
  mart.convert(wb.gene, "transcript", "public_name")(rownames(r.deconvolve))
q
g = intersect(rownames(embryo.timeseries), rownames(r.tophat))


# Utility to log-transform, and center (but not scale) a matrix.
f = function(x) {
  r = t( scale( log2(1 + t(x)), scale = FALSE))
  array(r, dim=dim(r), dimnames=dimnames(r))
}

#r = log2(1 + cbind(embryo.timeseries[g,], r.tophat[g,],
#  r.deconvolve[match(g, rownames(r.deconvolve)),]))
r = cbind(f(embryo.timeseries[g,]),
  f(r.tophat[g,]),
  f( r.deconvolve[ match(g, rownames(r.deconvolve)) , ] ))

r[is.na(r)] = 0
r = r[ apply(r, 1, var) > 0 , ]

system("mkdir -p git/unmix/seq/timing/embr_ts_comparison")

# write.table(round(r, 3), file="git/unmix/seq/timing/embr_ts_comparison/embr_ts.tsv",
#   sep="\t", row.names=TRUE, col.names=NA)

if (TRUE) {

  basefile = "git/unmix/seq/timing/embr_ts_comparison/embr_ts"
  nbproc = 8
  method = "correlation"
  link = "complete"
  hr <- hcluster(r, method = method, link = link, nbproc = nbproc)
  # XXX column clustering is arbitrary
  hc <- hcluster(t(r), method = method, link = link, nbproc = nbproc)
  hc$order = sort(hc$order)

  r2gtr(hr, file = paste(basefile, ".gtr", sep = ""))
#  r2atr(hc, file = paste(basefile, ".atr", sep = ""))

  r2cdt(hr, hc, r,
    file = paste(basefile, ".cdt", sep = ""))
}

# hclust2treeview(r, file="git/unmix/seq/timing/embr_ts_comparison/embr_ts.cdt")

