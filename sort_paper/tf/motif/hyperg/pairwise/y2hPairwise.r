# Runs the pairwise thing for the Y2H interactions.

source("git/data/name_convert.r")
source("git/sort_paper/tf/motif/hughes/motifInfo.r")
source("git/sort_paper/tf/motif/hyperg/pairwise/pairwiseMotif.r")

load("git/sort_paper/tf/motif/hughes/motifCluster.Rdata")

y2h = read.table("data/ppi/reeceHoyes2013_mmc7.tsv.gz",
  sep="\t", header=T, as.is=T)
y2h$bait.g = rename.gene.name.vector(y2h$bait.ORF.name)
y2h$prey.g = rename.gene.name.vector(y2h$prey.ORF.name)

# XXX this should be all motifs for which I ran FIMO,
# and it printed some output
all.scanned.motifs = unique(sub("_upstreamMotifCons.tsv.gz", "",
  c(list.files(paste0(motif.gene.dir, "/Ce_1.02")),
    list.files(paste0(motif.gene.dir, "/Dm_1.02")),
    list.files(paste0(motif.gene.dir, "/Mm_1.02")),
    list.files(paste0(motif.gene.dir, "/Hs_1.02")))))

# For a given organism, maps motifs to a canonical motif
# (e.g., what was used before.)
get.canonical.motifs = function(species) {
  a = cutree(hughes.motif.cluster[[species]], h = 0.01)
  a1 = a
  a1 = a1[ names(a1) %in% all.scanned.motifs ]
  canonical.motif = names(a1)[ match(a, a1) ]
  names(canonical.motif) = names(a)
  canonical.motif = canonical.motif[ !is.na(canonical.motif) ]
  canonical.motif
}

# XXX these should be the motifs I picked earlier
canonical.motif = hughes.motif.cluster[["Ce"]]$labels
names(canonical.motif) = canonical.motif
for(species in c("Dm", "Mm", "Hs")) {
  m = get.canonical.motifs(species)
  m = m[ !(names(m) %in% names(canonical.motif)) ]
  canonical.motif = c(canonical.motif, m)
}

# mapping from gene to motif, which for now is somewhat arbitrary
# motifs.1 = unique(c(nr.motifs, recursive=TRUE))
# motif.info.small = motif.info 

gene.to.motif =
  tapply(motif.info$motif.id, motif.info$gene, function(a) a[1])

y2h$bait.m = gene.to.motif[ y2h$bait.g ]
y2h$prey.m = gene.to.motif[ y2h$prey.g ]
y2h$bait.m.canonical = canonical.motif[ y2h$bait.m ]
y2h$prey.m.canonical = canonical.motif[ y2h$prey.m ]



motif.pairs = y2h[ , c("bait.m.canonical", "prey.m.canonical") ]
motif.pairs =
  motif.pairs[ (!is.na(motif.pairs[,1])) & (!is.na(motif.pairs[,2])) , ]
motif.pairs = unique(motif.pairs)

r = NULL
if (TRUE) {
  for(i in 1:nrow(motif.pairs)) {
    r = rbind(r,
      pairwise.motif.enrich(motif.pairs[i,1], motif.pairs[i,2]))
    if (i %% 25 == 0) {
      write.tsv(r,
        gzfile("git/sort_paper/tf/motif/hyperg/pairwise/y2hPairwise.tsv.gz"))
      cat(i, "\n")
    }
  }
  write.tsv(r,
    gzfile("git/sort_paper/tf/motif/hyperg/pairwise/y2hPairwise.tsv.gz"))
}

r$motif.dist.p.corr = p.adjust(r$motif.dist.p, method="fdr")

# annotate with what the original Y2H interactions were
colnames(r)[1] = "bait.m.canonical"
colnames(r)[2] = "prey.m.canonical"
r1 = r[ r$motif.dist.p.corr <= 0.1 , ]
y2h.pairwise.annotated = merge(y2h, r1)
write.tsv(y2h.pairwise.annotated,
  "git/sort_paper/tf/motif/hyperg/pairwise/y2hPairwiseAnnotated.tsv")

