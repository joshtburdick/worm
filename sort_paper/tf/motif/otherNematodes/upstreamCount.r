# Counts upstream motifs, for non-Ce species.

source("git/tf/motif/count/upstreamMotifCount.r")

# list of gene names to include (so that this is consistent across
# all the motifs)
gene.names = read.table(
    "git/data/seq/cb3/cb3_Ce_gene_upstream_3kb.bed",
    sep="\t", as.is=TRUE)[,4]

# various cutoffs to use in counting motifs
cutoffs = list(upstream.dist = c(1, 2, 3),
  conservation = c(0),
  score = c(30, 35, 40))

# Writes motif counts for one set of genes.
write.motif.counts = function(a) {
  cat(paste0("[writing motifs for ", a, "]\n"))
  motif.gene.dir =
    paste0("/media/jburdick/disk2/jburdick/distAndConservation_cb3/", a, "/")
  motifs = sub("_upstreamMotifCons.tsv.gz", "", list.files(motif.gene.dir))
  r = count.motifs.at.cutoffs(motifs,
    gene.names,
    cutoffs,
    get.motifs.from.file(motif.gene.dir),
    paste0("git/sort_paper/tf/motif/otherNematodes/upstreamCount/cb3/", a))
}

write.motif.counts("Ce_1.02")
write.motif.counts("Dm_1.02")
write.motif.counts("Mm_1.02")
write.motif.counts("Hs_1.02")