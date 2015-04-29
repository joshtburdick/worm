# Analysis of the qPCR results.

source("git/utils.r")

data.dir = "data/expression/qPCR/20150428 daf-19 hlh-6/"

r1 = read.tsv(paste0(data.dir, "results1.tsv"))
r2 = read.tsv(paste0(data.dir, "results2.tsv"))

# omit some controls
r2 = r2[ r2$"Sample Name" != "NORT" , ]
r2 = r2[ r2$"Sample Name" != "Water" , ]

# slight name tweaking

if (FALSE) {
r1$"Sample Name" = sub("daf1a", "daf19", r1$"Sample Name")
r1$"Sample Name" = sub("hih6", "hlh6", r1$"Sample Name")
r2$"Sample Name" = sub("hih6", "hlh6", r2$"Sample Name")
r1$"Target Name" = sub("hih6", "hlh6", r1$"Target Name")
r2$"Target Name" = sub("Phat5", "phat5", r2$"Target Name")
}

# Summarizes a table of RT results.
# Args:
#   r - a table of RT results
# Returns: table with columns
#   sample, replicated, target - description of the well
#   mean, sd, n - statistics for that combination
rt.summarize = function(r) {
  r = r[ ! (r$"Sample Name" %in% c("NORT", "Water")) , ]

  # XXX slight name tweaking
  r$"Sample Name" = sub("daf1a", "daf19", r$"Sample Name")
  r$"Sample Name" = sub("hih6", "hlh6", r$"Sample Name")
  r$"Target Name" = sub("Phat5", "phat5", r$"Target Name")
  ct = r$CT
  ct[ ct=="Undetermined" ] = NA
  ct = as.numeric(ct)

  # average technical replicates together
  r = aggregate(ct, list(sample.and.rep = r$"Sample Name", target = r$"Target Name"),
    function(x) {
      n = sum(!is.na(x))
      s = sd(x, na.rm=TRUE)
      c(mean=mean(x, na.rm=TRUE), sd=s, se = s / sqrt(n), n = n)
    },
    simplify=TRUE)
  r$ct = ifelse(r$x[,"n"] >= 3, r$x[,"mean"], NA)

  # split up sample name and replicate  
  re = "^(.+)(0\\d\\d\\d\\d\\d.?)$"
  m = regmatches(r$sample.and.rep, regexec(re, r$sample.and.rep))
  r$sample = sapply(m, function(a) a[2])
  r$rep = sapply(m, function(a) a[3])

  r = r[ !is.na(r$ct) , ]

  list(full = r, small = r[ , c("sample", "rep", "target", "ct") ])
}




s1 = rt.summarize(r1)$small
s2 = rt.summarize(r2)$small

