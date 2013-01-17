# Accuracy of sampling mean.

source("R/unmix/eval.r")

expr.cell = as.matrix(
  read.table("git/unmix/unmix_comp/data/exprCell.tsv",
  sep="\t", header=TRUE, row.names=1, as.is=TRUE))

load("git/unmix/comp_paper/sampling/plot/samplingStats.Rdata")

cat("cor. = ",
  mean(diag(cor(t(expr.cell), t(samplingStats[,"mean",])))), "\n")

cat("AUC = ",
  mean(auc(expr.cell >= 2000, samplingStats[,"mean",])), "\n")




