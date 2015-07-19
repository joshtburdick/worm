# Comparison of Y1H data with motif enrichment results.

source("git/utils.r")
source("git/data/name_convert.r")

# TF-cluster combined stats
if (FALSE) {
cluster.tf = read.tsv(gzfile("git/sort_paper/network/clusterTF.tsv.gz"))
cluster.tf = cluster.tf[ !is.na(cluster.tf[ , "TF-cluster corr." ]) , ]
cluster.tf = cluster.tf[ !is.na(cluster.tf[ , "Motif p" ]) , ]
cluster.tf = cluster.tf[ order(cluster.tf$"Motif p") , ]
cluster.tf = cluster.tf[ !duplicated(cluster.tf[ , c(1,2)]) , ]
}

# one of the clusterings
cl = {
  cl1 = read.table("git/cluster/hierarchical/hier.300.clusters/clusters.tsv",
    sep="\t", header=TRUE, as.is=TRUE)
  cl = cl1[,3]
  names(cl) = cl1[,2]
  cl
}

# Y1H data
y1h = read.delim(gzfile("data/tf/y1h/reeceHoyes2013_mmc3.tsv.gz"), as.is=T)

y1h$bait.name = rename.gene.name.vector(y1h$bait.ORF.name)
y1h$prey.name = rename.gene.name.vector(y1h$prey.ORF.name)
y1h$bait.cl = cl[ y1h$bait.name ]
y1h$prey.cl = cl[ y1h$prey.name ]
y1h = y1h[ (!is.na(y1h$bait.cl)) & (!is.na(y1h$prey.cl)) , ]

# simplified DBD names
simplified.dbd = {
  a = unique(c(y1h$bait.DBD, y1h$prey.DBD))
  names(a) = a
  for(tf.class in c("AT Hook", "bHLH", "HD", "HMG", "MYB",
      "Paired Domain", "WH", "ZF"))
    a[ grepl(paste0("^", tf.class), a) ] = tf.class
  a
}
y1h$bait.DBD = simplified.dbd[ y1h$bait.DBD ]
y1h$prey.DBD = simplified.dbd[ y1h$prey.DBD ]

# for each row of this, see if it has a Y1H entry
cluster.tf$has.y1h = paste(cluster.tf$TF, cluster.tf$Cluster) %in%
  paste(y1h$prey.name, y1h$bait.cl)

pdf("git/sort_paper/cluster/comparison/y1h.pdf")
par(mfrow=c(2,1))
hist(-log10(cluster.tf[ !cluster.tf$has.y1h , "Motif p" ]),
  xlim=c(0,10), breaks=500, col="grey",
  main="TF-gene pairs without Y1H interaction",
  xlab="-log10(p)")
hist(-log10(cluster.tf[ cluster.tf$has.y1h , "Motif p"]),
  xlim=c(0,10), breaks=500, col="grey",
  main="TF-gene pairs with Y1H interaction",
  xlab="-log10(p)")
wt = wilcox.test(-log10(cluster.tf[ ! cluster.tf$has.y1h , "Motif p" ]),
  -log10(cluster.tf[ cluster.tf$has.y1h , "Motif p"]),
  alternative="less")
print(wt)
print(wt$p.value)
dev.off()

# Number of interactions without, and with, a Y1H interaction
# (possibly a subset of Y1H interactions.)
num.y1h = function(cluster.tf, y1h, p.cutoff) {
  a = cluster.tf
  a$has.y1h = paste(cluster.tf$TF, cluster.tf$Cluster) %in%
    paste(y1h$prey.name, y1h$bait.cl)
  lp0 = -log10(a[!a$has.y1h,"Motif p"])
  lp1 = -log10(a[a$has.y1h,"Motif p"])
  c(no.y1h = mean(lp0 >= cutoff),
    with.y1h = mean(lp1 >= cutoff))
}

pdf("git/sort_paper/cluster/comparison/y1h_b.pdf",
  width=9, height=6)
par(mar=c(5,10,1,1)+0.1)
# par(mfrow=c(1,3))

prey.dbd = table(y1h$prey.DBD)
prey.dbd = prey.dbd[ prey.dbd >= 30 ]

for(cutoff in c(3)) {      # was c(2,3,4)) {

  r1 = num.y1h(cluster.tf, y1h, cutoff)
  names(r1) = c("No Y1H interaction", "All")

  r = NULL
  for(d in names(prey.dbd)) {
    print(d)
    r = c(r, num.y1h(cluster.tf, y1h[ y1h$prey.DBD == d,], cutoff)[2])
    names(r)[length(r)] = d
  }
  r[ is.na(r) ] = 0
  r = r[ order(r) ]

  r = c(r1["No Y1H interaction"], r, r1["All"])

  n = names(r)
  n[length(n)] = expression(bold("All"))

  barplot(r, beside=TRUE, horiz=TRUE, las=1,
    names.arg = n,
    col=c("white", rep("grey", length(r)-2), "#777777"),
#    main=paste("(cutoff =", cutoff, ")"),
    xlab=expression("Proportion of Y1H arcs with motif enrichment < 10"^{-3}))
#  legend("topleft",
#    legend=c("no Y1H", "with Y1H"),
#    fill=c("grey", "black"))
}

dev.off()

