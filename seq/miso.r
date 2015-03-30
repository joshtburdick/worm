# Utilities for dealing with MISO's output.

source("git/utils.r")

# Reads in summary files from one run.
read.miso = function(result.dir) {
  r = NULL
  for(f in list.files(paste0(result.dir, "/summary/"))) {
    as.type = sub(".miso_summary$", "", f)
    r1 = data.frame(as.type = as.type,
      read.table(paste0(result.dir, "/summary/", f),
        quote = "", header = TRUE, as.is = TRUE),
      stringsAsFactors = FALSE)
    r = rbind(r, r1)
  }

  r
}

# Reads all the files in a given directory
read.miso.dir = function(a) {
  r = NULL
  for(f in list.files(a)) {
    write.status(f)
    r1 = data.frame(file = f,
      read.miso(paste0(a, "/", f)),
      stringsAsFactors = FALSE)

    r = rbind(r, r1)
  }

  r
}



