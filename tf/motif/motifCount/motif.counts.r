# Computes motif counts.
# Currently "countMotifs.pl" does most of the work.

# the known motifs, 5kb upstream (without conservation filtering)
known.motifs.5kb  = function() {
m = read.table(gzfile("git/tf/motif/motifCount/motifs_5kbUpstream.tsv.gz"),
  sep="\t", header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE)

motif.list = read.table("data/tf/meme/motifList.tsv", sep="\t",
  header=TRUE, stringsAsFactors=FALSE)

# only using a subset of these
db.to.use = c("JASPAR_CORE_2009_insects",
  "JASPAR_CORE_2009_nematodes",
  "JASPAR_CORE_2009_vertebrates",
  "jolma2013")

motifs.to.use =
  motif.list[ motif.list$database %in% db.to.use, "id" ]
motifs.to.use = intersect(motifs.to.use, rownames(m))

t(m[ motifs.to.use, ])
}
# motif = known.motifs.5kb()

load("git/tf/motif/clusterUsingMotIV.Rdata")

# the motifs, in the intergenic (up to 5kb) regions
# XXX for whatever reason, the MX000045 motif is truncated.
# I have no idea why this is, but as it's part of PRODORIC
# (a prokaryotic database) anyway, for now I'm just ignoring the
# problem (it will be omitted anyway.) (The issue may be
# because that motif occurs a lot.)
known.motifs.1 = read.table(
  gzfile("git/tf/motif/motifCount/motifs_5kb_nogenes.tsv.gz"),
  sep="\t", header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE,
  fill=TRUE)
de.novo.motifs.1 = read.table(
  gzfile("git/tf/motif/motifCount/denovo_motifs_20130615_5kb_nogenes.tsv.gz"),
  sep="\t", header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE)

# restrict known motif set to those which were clustered
m1 = intersect(rownames(known.motifs.1), motif.clusters$label)

known.motifs = t(known.motifs.1[m1,])
de.novo.motifs = t(de.novo.motifs.1)

save(known.motifs, de.novo.motifs, file="git/tf/motif/motifCount/motif.counts.Rdata")

