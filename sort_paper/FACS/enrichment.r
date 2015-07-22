# The enrichments, pretty-printed.

source("git/utils.r")
source("git/sort_paper/plot/experimentRename.r")

r = read.tsv("git/cluster/readRatios.tsv")

colnames(r) = prettify.read.ratio.columns(colnames(r))

write.tsv(r, "git/sort_paper/FACS/enrichment.tsv")

