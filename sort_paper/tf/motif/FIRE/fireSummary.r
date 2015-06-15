# Summarizes results from FIRE.

source("git/utils.r")


# Summarizes the results in a directory.
# Args:
#   fire.result.dir - the directory of results
# Returns: a data frame
summarize.dir = function(fire.result.dir) {
  r = NULL

  for(d in list.files(fire.result.dir)) {
    cl = sub("_FIRE", "", d)
    write.status(cl)
    f = paste0(fire.result.dir, "/", d, "/DNA/", cl, ".summary")
    if (file.exists(f) && file.info(f)$size > 0) {

      # the summary info
      # XXX ignoring the last column
      a = read.table(f, sep="\t", header=FALSE, as.is=TRUE, fill=TRUE)[,c(1:12)]
      colnames(a) = c("optimizedMotif", "x1", "x2", "mi", "x3", "zScore",
        "robustness", "x4", "seed", "x5", "x6", "conservationIndex")
#      a$cluster = cl
      r = rbind(r, cbind(cluster=cl, a))
    }
  }

  r[,1] = as.character(r[,1])
  r
}


facs.summary = summarize.dir("git/sort_paper/tf/motif/FIRE/output/facs/")
cat(paste("\nFACS: ", nrow(facs.summary), " motifs found near ",
  length(unique(facs.summary$cluster)), " clusters\n"))
spencer.summary = summarize.dir("git/sort_paper/tf/motif/FIRE/output/spencerEmbryonic/")
cat(paste("\nSpencer expression: ", nrow(spencer.summary), " motifs found near ",
  length(unique(spencer.summary$cluster)), " clusters\n"))

output.dir = "git/sort_paper/tf/motif/FIRE/fireSummary/"
system(paste("mkdir -p ", output.dir))
write.tsv(facs.summary, paste0(output.dir, "facs.tsv"))
write.tsv(spencer.summary, paste0(output.dir, "spencer.tsv"))


