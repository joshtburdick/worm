# Purity of sorting, particularly comparing
# ungated to singlets.

r = as.matrix(read.table("R/unmix/sort_paper/seq/quant/readsPerMillion.tsv",
  header=TRUE, row.names=1, as.is=TRUE))

lr = log2(1+r)

# Differences between unsorted and singlets.
ungated.singlet.diff =
  ((lr[,"pha.4_ungated"] + lr[,"cnd.1_ungated"])
  - (lr[,"pha.4_singlets"] + lr[,"cnd.1_ungated"]))

pdf("git/unmix/seq/FACS/ungated_and_singlets.pdf", width=7, height=7)
par(mfrow=c(2,2))
hist(lr[,"cnd.1_ungated"] - lr[,"cnd.1_singlets"],
  main="cnd-1 ungated - cnd-1 singlets",
  xlab="Difference in log-transformed coverage",
  breaks=50, col="grey")
hist(lr[,"pha.4_ungated"] - lr[,"pha.4_singlets"],
  main="pha-4 ungated - pha-4 singlets",
  xlab="Difference in log-transformed coverage",
  breaks=50, col="grey")
smoothScatter(lr[,"cnd.1_ungated"] - lr[,"cnd.1_singlets"],
  lr[,"pha.4_ungated"] - lr[,"pha.4_singlets"],
  xlab="cnd-1 ungated - cnd-1 singlets",
  ylab="pha-4 ungated - pha-4 singlets")

cat("correlation of ungated - singlets across both = ",
  cor(lr[,"cnd.1_ungated"] - lr[,"cnd.1_singlets"],
    lr[,"pha.4_ungated"] - lr[,"pha.4_singlets"]),
  "\n")

hist(ungated.singlet.diff,
  main="Combined ungated - singlets",
  xlab="Difference in log-transformed coverage",
  breaks=50, col="grey")

dev.off()

write.table(names(ungated.singlet.diff)[ungated.singlet.diff > 2],
  file="git/unmix/seq/FACS/ungated_high.txt",
  row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(names(ungated.singlet.diff)[ungated.singlet.diff < -1.5],
  file="git/unmix/seq/FACS/ungated_low.txt",
  row.names=FALSE, col.names=FALSE, quote=FALSE)



# Does many t-tests between x and y. Not radically fast.
# Returns: data.frame containing
#   t - the t-value
#   p - the uncorrected p-value
t.test.many = function(x, y) {
  r = NULL

  for(i in 1:nrow(x)) {
    if (i && 1000 == 0) {
      cat(i, "")
    }
    if (var(c(x[i,], y[i,])) < 1e-2) {
      r = rbind(r, c(t=NA, p=NA))
    }
    else {
      a = t.test(x[i,], y[i,])
      r = rbind(r, c(a$statistic[1], p = a$p.value))
    }
  }

  rownames(r) = rownames(x)
  r
}




ungated.singlets.t = data.frame(t.test.many(lr[,c(18,24)], lr[,-c(18,24)]))
ungated.singlets.t[ is.na(ungated.singlets.t$p), "p" ] = 1
ungated.singlets.t$adjusted.p = p.adjust(ungated.singlets.t$p, method="hochberg")
ungated.singlets.t =
  ungated.singlets.t[order(ungated.singlets.t$adjusted.p, decreasing=FALSE),]
write.table(ungated.singlets.t, file="git/unmix/seq/FACS/ungated_singlets_t.tsv",
  sep="\t", row.names=TRUE, col.names=NA)

# the p-values here don't seem that significant
if (FALSE) {
write.table(rownames(ungated.singlets.t[!is.na(ungated.singlets.t$adjusted.p) &
    ungated.singlets.t$adjusted.p <= 0.05 &
    ungated.singlets.t$t > 0,]),
  file="git/unmix/seq/FACS/ungated_high.txt",
  row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(rownames(ungated.singlets.t[!is.na(ungated.singlets.t$adjusted.p) &
    ungated.singlets.t$adjusted.p <= 0.05 &
    ungated.singlets.t$t < 0,]),
  file="git/unmix/seq/FACS/ungated_low.txt",
  row.names=FALSE, col.names=FALSE, quote=FALSE)
}

