# Plots mean and s.d. of when genes are expressed.

source("git/utils.r")

gene.expr.time =
  read.tsv("git/sort_paper/FACS/timing/geneExprTime.tsv")

sd.cutoff = quantile(gene.expr.time$sd, p=0.25)
cat("0.25 quantile of s.d. =", sd.cutoff, "\n")

pdf("git/sort_paper/FACS/timing/geneExprMeanAndSD.pdf",
  width=6, height=6)
smoothScatter(gene.expr.time$mean, gene.expr.time$sd,
  main="Time when genes are expressed",
  xlab="mean of time expressed (minutes)",
  ylab="s.d. of time expressed (minutes)")
abline(h=sd.cutoff, col = "#00000080")
dev.off()

