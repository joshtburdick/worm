# Counts of TFs associated with particular clusters.

source("git/utils.r")

ctf = read.tsv(gzfile("git/sort_paper/network/clusterTF.tsv.gz"))

a = ctf[ctf$TF=="pha-4" & ctf$"Motif enrich" >= 2,]
a = a[ , c("Cluster", "TF", "TF-cluster corr.", "Motif", "Motif p")]
a = a[ order(a$"Motif p") , ]
a = a[ !duplicated(a[,c("Cluster", "TF")]) , ]

# picking a cutoff
a = a[ a$"TF-cluster corr." >= 0.4 & a$"Motif p" <= 1e-4 , ]
print(a)

# ChIP
cl = a$Cluster
pha4.chip = ctf[ctf$Cluster %in% cl , ]
pha4.chip = pha4.chip[ order(pha4.chip$"ChIP p") , ]
pha4.chip = pha4.chip[ !duplicated(pha4.chip$Cluster) , ]

print(pha4.chip)
