
source("git/utils.r")
source("git/plot/seqLogoSVG.r")

load("git/tf/motif/meme.format.pwm.Rdata")

# Writes one motif as an SVG file.
# Args:
#   pwm - a position-weight matrix
#   output.file - SVG file to write
# Side effects: writes a file
write.motif.svg.orig = function(pwm, output.file) {
  svg(output.file, width=60, height=10)
  par(mar=c(0,0,0,0))
  seqLogo(pwm, xaxis=FALSE, yaxis=FALSE, xfontsize=0, yfontsize=0)
 
  dev.off()
}

hughes.motif.table = read.table(
  "data/tf/hughes/TF_Information_all_motifs.tsv",
  sep="\t", header=TRUE, as.is=TRUE)

hughes.motifs = unique(hughes.motif.table$Motif_ID)

# hughes.motifs = c("M6986_1.00", "M1112_1.00")

output.dir = "git/tf/motif/plot/motifSvg_hughes/"
system(paste("mkdir -p", output.dir))

# for (m in "GATA3_DBD") {
for (m in intersect(hughes.motifs, names(meme.format.pwm))) {

  write.status(m)
  write.motif.svgz(meme.format.pwm[[m]],
    paste0(output.dir, "/", m, ".svg"))
}

cat("\n")

