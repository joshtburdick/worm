# The deconvolved embryo timeseries, with some renaming.

source("git/utils.r")
source("git/data/biomart_utils.r")

# deconvolved data from modENCODE
r.deconvolve = as.matrix(read.table(gzfile("data/modencode/N2_EE_50_mod1.csv.gz"),
  sep=",", header=FALSE, row.names=1))
colnames(r.deconvolve) = c("t.000", "t.060", "t.120", "t.150", "t.180",
  "t.240", "t.330", "t.390", "t.420", "t.480",
  "t.540", "t.570", "t.600", "t.630", "t.660")

# change row names to "public_name"
rownames(r.deconvolve) =
  mart.convert(wb.gene, "transcript", "public_name")(rownames(r.deconvolve))

# normalize to "reads per million" (not subtracting out ribosomal RNA,
# as it's not included in these counts)
r.deconvolve = t( t(r.deconvolve) / (apply(r.deconvolve, 2, sum) / 1e6) )

write.tsv(r.deconvolve, "git/unmix/seq/timing/deconvolved_embryo_ts.tsv")

