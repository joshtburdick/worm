
library("flashClust")
library("amap")

# get data (FIXME should be in a separate function)
count.path = "git/unmix/seq/quant/readsPerMillion"

r1 = as.matrix(
  read.table(paste(count.path, "readsPerMillion_Murray_050912.tsv", sep="/"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

r2 = as.matrix(
  read.table(paste(count.path, "readsPerMillion_092812.tsv", sep="/"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

r = log2(1 + cbind(r1, r2))
r = r[ , order(colnames(r)) ]

# remove constant rows
r = r[ apply(r, 1, var) > 0 , ]

# if (FALSE) {
print(date())
  a = cor(t(r))
print(date())
  a = 1 - a
  gc()
print(date())
  # XXX for whatever inane reason, this takes like eight minutes.
  # not worrying about this for now.
  d = as.dist(a)
print(date())
  a = NULL
  gc()
save(d, file="/var/tmp/exprCor.Rdata")
print(date())
# }

# h = hclust(d, method="complete")
# print(date())
save(h, file="/var/tmp/exprClustCorComplete.Rdata")

