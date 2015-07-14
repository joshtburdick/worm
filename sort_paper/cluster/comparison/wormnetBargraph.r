# Bar graph of Wormnet enrichments.

source("git/utils.r")

w = read.tsv("git/sort_paper/cluster/comparison/wormnet.tsv")
w = w[ order(
  rownames(w)=="combined",
  w$times.more.than.random) , ]
w["combined","description"] = expression(bold("All"))

pdf("git/sort_paper/cluster/comparison/wormnetBargraph.pdf",
  width=11, height=6)
par(mar=c(5,29,0.7,0.7)+0.1)

barplot(w$times.more.than.random, horiz=TRUE,
  names.arg = w$description, las=1,
    col=ifelse(rownames(w)=="combined", "#777777", "grey"),
  xlab="Fold enrichment over random")

dev.off()

