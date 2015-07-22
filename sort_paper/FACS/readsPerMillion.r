# The RPM counts, pretty-printed.

source("git/utils.r")
source("git/sort_paper/plot/experimentRename.r")

r = read.tsv("git/cluster/readsPerMillion.tsv")

r = r[ , !grepl("^(HS|RNAi|lowInput|normalInput)", colnames(r)) ]

colnames(r) = prettify.read.ratio.columns(colnames(r))

write.tsv(r, "git/sort_paper/FACS/readsPerMillion.tsv")

