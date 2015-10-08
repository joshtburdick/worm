# Make-like utilities for caching results of scripts.

source("git/utils.r")

# Sources an R script, in a separate process.
# (This is intended to avoid storing variables computed when
# sourcing previous files; admittedly it's slower.)
# Args:
#   source.file - the script to run
#   save.env - if TRUE, save environment in an .Rdata file.
# Side effects: whatever the script does (and
#   possibly writes an .Rdata file.)
source.separate.process = function(source.file, save.env) {
  rdata.file = sub("\\.[Rr]$", ".Rdata", source.file)

  # construct a command to source the file, and maybe save results
  r = paste0('source("', source.file, '"); ')
  if (save.env)
    r = paste0(r,
      'save.image(file="', rdata.file, '"); ')

  cmd = paste0("R --quiet --no-save -e '", r, "' > /dev/null")

  print(cmd)
}

# Version of source(), with caching.
src = function(f) {



  cacheFile = sub("\\.[Rr]$", ".Rdata", f1)
#  if (



}

# Ensures that a given file is present and up-to-date.
#   If f ends in .r or .R, it's assumed to be an R script;
# use() will check for a corresponding up-to-date .Rdata file.
# Otherwise, use() will look for a script with the same name,
# and assume that script generated the output file.
# (by changing the extension to .r). In either case, if the
# script (or any dependency) is newer than the corresponding
# source file, the script will be re-run.
# This won't create an .Rdata file .
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

  
  # guess the script name (for now, only searching for ".R";
  # ".r" may be more standard)
  s = sub("\\.[^\\.]+$", ".[Rr]", f1)


  if (is.script) {
  
  }
  else {


  }

  f
}

# Gets the dependencies of a given R source file.
# Note that this is defined _very_ naively, by simply
# scanning for strings of the form
#
#   src("...filename...")
#
# (or use() or source()), all on one line.
#
# Args:
#   source.file - name of the file to scan
# Returns: list of dependent files.
# XXX not yet working
get.dependencies.one.file = function(source.file) {

  # read all lines in file
  s = scan(source.file, sep="\n", what=list(""), quiet=TRUE)[[1]]

  f = sapply(z, function(a)
    regexp.capture(a, "(?:source|src|use)\\(\"([^\"]+)\"\\)"))

  f

}

# Recursively gets the dependencies.
rebuild = function(source.file) {


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

