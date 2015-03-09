# Summarizes results of shuffling motifs.

source("git/utils.r")

result.dir = "git/tf/motif/shuffle/shuffledCompare/"

motifs = sub(".Rdata", "", list.files(result.dir))


# Gets summary stats for one motif.
get.summary = function(m) {
  shuffled.motif.dist = NULL
  load(paste0(result.dir, "/", m, ".Rdata"))
  c(mean(shuffled.motif.dist), quantile(shuffled.motif.dist))
}

r = matrix(NA, nrow=length(motifs), ncol=6)
rownames(r) = motifs
colnames(r) = c("mean", "q.0", "q.25", "q.50", "q.75", "q.100")

for(m in motifs) {
  write.status(m)
  r[m,] = get.summary(m)
}

write.tsv(r, "git/tf/motif/shuffle/shuffledCompareStats.tsv")

