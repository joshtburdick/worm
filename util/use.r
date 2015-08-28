# Utilities for incrementally computing things.
# For now, doesn't automatically detect stale dependencies.

# Ensures that a given file is present and up-to-date.
#   If f ends in .r or .R, it's assumed to be an R script;
# use() will check for a corresponding up-to-date .Rdata file.
# Otherwise, use() will look for a script with the same name,
# and assume that script generated the output file.
# (by changing the extension to .r). In either case, if the
# script (or any dependency) is newer than the corresponding
# source file, the script will be re-run.
# Args:
#   f - the name of a file which is needed.
#     use() behaves differently depending on what f is.
# Returns: f. This may be discarded, or can be used,
# e.g. in something like:
#   z = read.table(use("somefile.tsv"), sep="\t")
# XXX that's not yet implemented.
# Side effects: if f is a script, reads in its saved definitions.
use = function(f) {
  # first, remove a filename suffix which "looks like" compression
  f1 = sub("\\.(bz|bz2|gz|lz|lzma|rz|sz|xz|z|Z)$", "", f)

  


  # guess the script name (for now, only searching for ".r";
  # ".R" may be more standard)
  s = sub("\\.[^\\.]+$", ".r", f1)




  if (is.script) {
  
  }
  else {


  }

  f
}

# Gets the dependencies of a given R source file.
# Note that this is defined _very_ naively, by simply
# scanning for strings (all on one line) of the form
#
#   use("...filename...")
#
# Args:
#   source.file - name of the file to scan
# Returns: list of dependent files.
# XXX not yet working
get.dependencies = function(source.file) {
  # read all lines in file
  s = scan(source.file, sep="\n", what=list(""), quiet=TRUE)[[1]]

  f = sapply(z, function(a) regexp.capture(a, "use\\(\"([^\"]+)\"\\)"))


  f

}

# Utility which gets the name of the current script, with
# .r or .R replaced with something else.
# Useful for writing an output file named similarly to a script.
# Args:
#   extension - the extension to add on to the filename
#     (which should include the ".").
# Returns: name of the current script, but with the
#   suffix replaced.
same.name = function(extension) {
# from
# http://stackoverflow.com/questions/13672720/...
# r-command-for-setting-working-directory-to-source-file-location
  script.name = sys.frame(1)$ofile
  base.name = sub("\\.[rR]$", "", script.name)
  paste0(base.name, extension)
}

