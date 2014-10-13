# Reads in all motifs which are in MEME format.

library("MotIV")

# hack to read a MEME format file directly.
read.meme.file = function(f) {
  cat("reading", f, "\n")

  # XXX should use a proper temp file, and, like, check for errors

#  system(paste("/home/jburdick/meme/bin/meme2meme ", f, " > /var/tmp/meme.txt"))
#  system(paste("/home/jburdick/gcb/git/perl/motif/meme2TRANSFAC.pl /var/tmp/meme.txt /var/tmp/meme.transfac.txt")
  system(paste("/home/jburdick/gcb/git/perl/motif/meme2TRANSFAC.pl ", f, " /var/tmp/meme.transfac.txt"))
  r = readPWMfile("/var/tmp/meme.transfac.txt")
  r
}

# Reads in all MEME files in a directory.
read.meme.files.in.dir = function(meme.dir) {
  r = NULL
  for(f in list.files(meme.dir)) {
    f1 = paste(meme.dir, "/", f, sep="")
    r = c(r, read.meme.file(f1))
  }
  r
}


meme.format.pwm = {

  known.motif.meme.file = paste("data/tf/meme/motif_databases/",
    c("jolma2013.meme",
    "JASPAR_CORE_2009_insects.meme",
    "JASPAR_CORE_2009_nematodes.meme",
    "JASPAR_CORE_2009_vertebrates.meme"), sep="")

#  r = read.meme.files.in.dir("git/tf/motif/meme_file/")

  # XXX this is more complicated than presumably it needs to be
  r = c(read.meme.files.in.dir("git/tf/motif/meme_file/"),
    read.meme.file(known.motif.meme.file[[1]]),
    read.meme.file(known.motif.meme.file[[2]]),
    read.meme.file(known.motif.meme.file[[3]]),
    read.meme.file(known.motif.meme.file[[4]]),
    read.meme.file("git/tf/motif/hughes_motif.meme"))

  number.columns = function(x) {
# print(dim(x))
    colnames(x) = c(1:ncol(x))
    x
  }

  lapply(r, number.columns)
}

save(meme.format.pwm, file="git/tf/motif/meme.format.pwm.Rdata")

