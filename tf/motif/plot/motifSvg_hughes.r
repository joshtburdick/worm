
source("git/utils.r")
source("git/plot/seqLogoSVG.r")

load("git/sort_paper/tf/motif/hughes/motifMatrix.Rdata")

output.dir = "git/tf/motif/plot/motifSvg_hughes/"
system(paste("mkdir -p", output.dir))

# Writes a set of motifs.
write.motifs = function(pwm) {
  for (m in names(pwm)) {
    write.status(m)
    write.motif.svgz(pwm[[m]], paste0(output.dir, "/", m, ".svg"))
  }
  cat("\n")
}

cat("Ce\n")
write.motifs(hughes.motif.matrix[["Ce"]])

for(org in c("Dm", "Hs", "Mm")) {
  cat(org, "\n")
  m = hughes.motif.matrix[[org]]
  motif.list = read.table(
    paste0("git/sort_paper/tf/motif/hughes/", org, "_motif_nr.txt"),
      as.is=TRUE)[,1]
  write.motifs(hughes.motif.matrix[[ org ]][ motif.list ])
}

