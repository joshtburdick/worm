# Summarizes MISO results.

source("git/seq/miso.r")
source("git/utils.r")

r = read.miso.dir("/home/jburdick/tmp/miso")

write.tsv(r, file=gzfile("git/sort_paper/splicing/MISO/embryo_ts/misoSummary.tsv.gz"))

