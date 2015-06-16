# Statistics about the table of clusters and TFs.
# (And some checks on individual clusters.)

source("git/utils.r")

ctf = read.tsv(gzfile("git/sort_paper/network/clusterTF.tsv.gz"))

# very basic positive controls
ctf1 = ctf[ ctf$Cluster %in% c("30", "52", "286") &
  ctf$TF %in% c("pha-4", "let-381", "daf-19") , ]
ctf1 = ctf1[ ! duplicated(ctf1[,c("Cluster", "TF")]) , ]
print(ctf1)

# version of this with NAs removed (and fewer columns)
ctf2 = ctf[ , c("Cluster", "TF", "TF-cluster corr.", "Motif p", "ChIP p") ]
ctf2$"Motif p"[ is.na(ctf2$"Motif p") ] = 1
ctf2$"ChIP p"[ is.na(ctf2$"ChIP p") ] = 1

# table of number
num.clusters.with.tfs = NULL
for(cor.cutoff in c(-1, -0.5, -0.4, 0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
  for(motif.p.cutoff in c(1, 0.05, 1e-3, 1e-5, 1e-10))
    for(chip.p.cutoff in c(1, 0.05, 1e-3, 1e-5, 1e-10)) {
      write.status(paste(cor.cutoff, motif.p.cutoff, chip.p.cutoff))
      a = ctf2[ ctf2$"TF-cluster corr." >= cor.cutoff &
        ctf2$"Motif p" <= motif.p.cutoff &
        ctf2$"ChIP p" <= chip.p.cutoff , c("Cluster", "TF") ]
      a = unique(a)
      b = c("Correlation cutoff" = cor.cutoff,
        "Motif p cutof" = motif.p.cutoff,
        "ChIP p cutoff" = chip.p.cutoff,
        num.clusters = length(unique(a$Cluster)),
        num.tfs = length(unique(a$TF)),
        num.pairs = nrow(a))
      num.clusters.with.tfs = rbind(num.clusters.with.tfs, b)
    }

write.tsv(num.clusters.with.tfs, "git/sort_paper/network/numClustersWithTFs.tsv")

