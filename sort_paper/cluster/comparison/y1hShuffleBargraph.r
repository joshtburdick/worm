# Bar graph of shuffled Y1H enrichments.

source("git/utils.r")
source("git/plot/utils.r")

w = read.tsv("git/sort_paper/cluster/comparison/y1h_graph_shuffle_2_1e5.tsv")
w = w[ order(
  rownames(w)=="All",
  w$overlap.enrich) , ]
w$description = rownames(w)
w["All","description"] = expression(bold("All"))

# XXX omitting "AT hook" and "HMG", because they're so sparse
w = w[ !(rownames(w) %in% c("AT Hook", "HMG")) , ]

pdf("git/sort_paper/cluster/comparison/y1hShuffleBargraph2.pdf",
  width=6, height=4)
par(mar=c(5,5,2,6)+0.1)

b = barplot(w$overlap.enrich, horiz=TRUE,
  names.arg = w$description, las=1,
    col=ifelse(rownames(w)=="All", "#777777", "grey"),
  xlab="Fold enrichment over random")

x = max(w$overlap.enrich) + 0.5
for(i in 1:nrow(w)) {
  text(x, b[i], xpd=NA, adj=0,
    labels=expr.format(expression("p " * p),
      list(p=format.p(w$p.corr[i], include.equals.sign=TRUE))))
}

dev.off()

