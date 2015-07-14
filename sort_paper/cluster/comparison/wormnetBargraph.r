# Bar graph of Wormnet enrichments.

source("git/utils.r")

w = read.tsv("git/sort_paper/cluster/comparison/wormnet.tsv")
w = w[ order(w$times.more.than.random) , ]
w["combined","description"] = "combined"

pdf("git/sort_paper/cluster/comparison/wormnetBargraph.pdf",
  width=8, height=11)
par(mar=c(30,4,4,1)+0.1)

barplot(w$times.more.than.random,
  names.arg = w$description, las=2,
    col=ifelse(w$description=="combined", "#777777", "grey"),
  ylab="Fold enrichment over random")

dev.off()

