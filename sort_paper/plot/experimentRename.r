# Renames some columns (ideally this would happen
# earlier, but that might break some analyses.)

source("git/utils.r")

# experimentNames = read.tsv("git/unmix/seq/quant/experimentNames.tsv")

# Renames the read ratio columns.
# Args: names of columns
# Returns: names of columnes, renamed.
prettify.read.ratio.columns = function(cols) {
  r = cols

  r = sub("8/19", "rep. 1", r)
  r = sub("12/14", "rep. 2", r)
  r = sub("1/4", "rep. 3", r)

  r = sub("5/9", "rep. 1", r)
  r = sub("9/1", "rep. 2", r)
  r = sub("12/9", "rep. 3", r)

  i = grepl("^t\\.", cols)
  r[i] = as.character(as.numeric(sub("t\\.", "", r[i])))
  r
}



