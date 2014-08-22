# Creates SVG files of each motif.
# FIXME: possibly make these smaller?

library("seqLogo")
library("MotIV")

source("git/utils.r")

load("git/tf/motif/meme.format.pwm.Rdata")

# running this on the filtered (non-redundant) motifs
motif.gene.dir = "git/cluster/motif/distAndConservation/5kb/"
known.motifs = {
  r = list.files(motif.gene.dir)
  sub("_upstreamMotifCons.tsv.gz", "", r)
}
motif.filter = read.tsv("git/tf/motif/motifFilter.tsv")
known.motifs.small =
  intersect(known.motifs, motif.filter$canonical.name)

# Writes one motif as an SVG file.
# Args:
#   pwm - a position-weight matrix
#   output.file - SVG file to write
# Side effects: writes a file
write.motif.svg = function(pwm, output.file) {
  svg(output.file, width=60, height=10)
  par(mar=c(0,0,0,0))
  seqLogo(pwm, xaxis=FALSE, yaxis=FALSE, xfontsize=0, yfontsize=0)
 
  dev.off()
}

output.dir = "git/tf/motif/plot/motifSvg/"
system(paste("mkdir -p", output.dir))

for (m in intersect(rownames(motif.filter), known.motifs)) {   #  known.motifs.small) {
  write.status(m)
  write.motif.svg(meme.format.pwm[[m]],
    paste0(output.dir, "/", m, ".svg"))
}

