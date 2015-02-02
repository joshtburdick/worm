# Imports expression cluster data.

source("git/utils.r")

r = read.tsv(gzfile(
  "data/wormbase/ExpressionCluster/ExprClusterTable_WS244.tsv.gz"))
expr.cluster.descr = r[,1]
names(expr.cluster.descr) = rownames(r)

# in cases where clusters all have the same description,
# tack on the cluster ID
descr.count = table(expr.cluster.descr)
non.unique.descr = names(descr.count)[ descr.count > 1 ]
i =  expr.cluster.descr %in% non.unique.descr
expr.cluster.descr[ i ] = paste(expr.cluster.descr[i],
  names(expr.cluster.descr)[i])

r = read.tsv(gzfile(
  "data/wormbase/ExpressionCluster/GeneExprCluster_WS244.tsv.gz"))
# r = r[1:5,]  # XXX testing

# Utility converting a named list of vectors to
# a table.
list.to.table = function(a) {
  r = NULL

  # XXX this is not terrifically efficient
  for(i in 1:length(a)) {
    write.status(i)
    name = names(a)[i]
    r = rbind(r, cbind(name = name, x = a[[i]]))
  }

  data.frame(r, stringsAsFactors=FALSE)
}

a = strsplit(r[,2], ",")
names(a) = r[,1]
expr.cluster = list.to.table(a)
colnames(expr.cluster) = c("gene", "group")


save(expr.cluster, expr.cluster.descr,
  file="git/data/wormbase/expr.cluster.Rdata")

