# Comparison of how much our sort fractions overlap,
# compared with e.g. Spencer's.

source("git/utils.r")
source("git/plot/utils.r")
source("git/plot/plot_expr.r")

# the sort matrix (unnormalized)
source("git/sort_paper/unmix/sortMatrix.r")

# just show the leaf nodes
m.leaf = m.unnormalized %*% cell.to.leaf.matrix

# XXX
colnames(m.leaf)[colnames(m.leaf)=="NA."] = "P0"

# omitting double-negatives, for consistency with other rows
m.leaf = m.leaf[ rownames(m.leaf) != "ceh-6 (-) hlh-16 (-)" , ]

spencer.m = as.matrix(read.tsv(
  "git/sort_paper/validation/spencerSortMatrix.tsv"))
spencer.m.leaf = spencer.m %*% cell.to.leaf.matrix

# quantize these
facs.on.off = 1 * (m.leaf > 0)
spencer.on.off = 1 * (spencer.m.leaf > 0)

# overlap between all pairs of markers
overlap = function(m) {
  a = m %*% t(m)
  diag(a) = NA
  a = as.vector(a)
  a[ !is.na(a) ]
}

cat(paste0("FACS num cells    mean = ",
  mean(apply(facs.on.off,1,sum)),
  "    median = ", median(apply(facs.on.off,1,sum)), "\n"))
cat(paste0("Spencer num cells    mean = ",
  mean(apply(spencer.on.off,1,sum)),
  "    median = ", median(apply(spencer.on.off,1,sum)), "\n"))

cat(paste0("FACS num. fractions    mean = ",
  mean(apply(facs.on.off,2,sum)),
  "    median = ", median(apply(facs.on.off,2,sum)), "\n"))
cat(paste0("Spencer num. fractions   mean = ",
  mean(apply(spencer.on.off,2,sum)),
  "    median = ", median(apply(spencer.on.off,2,sum)), "\n"))

facs.overlap = overlap(facs.on.off)
spencer.overlap = overlap(spencer.on.off)

cat(paste0("FACS mean overlap = ", mean(facs.overlap), "\n"))
cat(paste0("Spencer mean overlap = ", mean(spencer.overlap), "\n"))

wt = wilcox.test(facs.overlap, spencer.overlap, alternative="greater")


pdf("git/sort_paper/unmix/sortMatrix/spencerOverlapCompare.pdf",
  width=4.5, height=6)
par(mfrow=c(2,1))
b = seq(0, 110, 2)
hist(facs.overlap, breaks=b, col="#ff0000a0",
  main="Pairs of markers used for flow sorting",
  xlab="Number of cells in overlap")
hist(spencer.overlap, breaks=b, col="#ffff00a0",
  main="Pairs of embryonic tissues\nfrom (Spencer et al. 2011)",
  xlab="Number of cells in overlap")
# XXX
x = expression("Wilcoxon p" * 1)
x[[1]][[3]] = format.p(wt$p.value)[[1]]
mtext(x, line=-1, cex=0.9)
dev.off()

