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
motif = known.motifs.5kb()


# the motifs, in the intergenic (up to 5kb) regions
#known.motifs.1 = read.table(
#  gzfile("git/tf/motif/motifCount/motifs_5kb_nogenes.tsv.gz"),
#  "git/tf/motif/motifCount/foo.tsv",
#  sep="\t", header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE)
#de.novo.motifs.1 = read.table(
#  gzfile("git/tf/motif/motifCount/denovo_motifs_20130615_5kb_nogenes.tsv.gz"),
#  sep="\t", header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE)

# known.motifs = 


