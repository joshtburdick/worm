# Reads in a MEME-format file of PWMs.

library("MotIV")

# hack to read a MEME format file directly.
read.meme.file = function(f) {
  cat("reading", f, "\n")

  # XXX error checking is a little hokey
  p = pipe(paste("git/perl/motif/meme2TRANSFAC.pl ", f), open="r")
  r = readPWMfile(p)
  close(p)

  number.columns = function(x) {
    colnames(x) = c(1:ncol(x))
    x
  }

  lapply(r, number.columns)
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

