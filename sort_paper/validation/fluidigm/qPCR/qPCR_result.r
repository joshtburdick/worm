# Analysis of the qPCR results.

source("git/utils.r")

data.dir = "data/expression/qPCR/20150428 daf-19 hlh-6/"

r1 = read.tsv(paste0(data.dir, "results1.tsv"))
r2 = read.tsv(paste0(data.dir, "results2.tsv"))

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
  r$"Target Name" = sub("hih6", "hlh6", r$"Target Name")
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

  # omit cases in which the CT wasn't defined
  r = r[ !is.na(r$ct) , ]

  list(full = r, small = r[ , c("sample", "rep", "target", "ct") ])
}

# Matches up a particular control name in all the experiments.
matched.control = function(s, control.name ) {
  s.control = s[ s$target==control.name , ]
  rownames(s.control) = paste(s.control$sample, s.control$rep)
  s.control[ paste(s$sample, s$rep) , "ct" ]
}

# Compares

# Args:
#   s - a data.table of qPCR results
#   sample, target - the sample and target to get numbers for
# Returns: vector of ddct stats, and (unadjusted) p-value
rt.compare.1 = function(s, sample, target, dct.column) {
  dct.sample = s[ s$sample==sample & s$target==target , dct.column ]
  dct.N2 = s[ s$sample=="N2" & s$target==target , dct.column ]
  ddct = mean(dct.sample, na.rm=TRUE) - mean(dct.N2, na.rm=TRUE)
  r = t.test(dct.sample, dct.N2, alternative="two.sided", var.equal=FALSE)
  data.frame(sample = sample, target = target, ddct = ddct,
    t = r$statistic, p = r$p.value,
    n.sample=length(dct.sample), n.N2=length(dct.N2), stringsAsFactors=FALSE)
}



rt.compare = function(s, dct.column) {
  r = NULL

  a = unique(s[s$sample != "N2", c("sample", "target")])
  a = a[ order(a$sample, a$target), ]

  for(i in 1:nrow(a)) {
    r = rbind(r, rt.compare.1(s, a[i,"sample"], a[i,"target"], dct.column))
  }
  rownames(r) = NULL

  r
}


s1 = rt.summarize(r1)$small

# remove an outlier (even though it nominally worked)
s1 = s1[ !(s1$sample=="N2" & s1$rep=="022415c" & s1$target=="ama1") , ]
s2 = rt.summarize(r2)$small
# slight name adjustment
s2$rep = sub("022414c", "022415c", s2$rep)

s = rbind(s1, s2)

if (TRUE) {
s$ct.ama1 = matched.control(s, "ama1")
s$ct.arf3 = matched.control(s, "arf3")
s$dct.ama1 = s$ct.ama1 - s$ct
s$dct.ama1arf3 = (s$ct.ama1 + s$ct.arf3) / 2 - s$ct


write.tsv(s, "git/sort_paper/validation/fluidigm/qPCR/qPCR_result.tsv")

print(rt.compare(s, "dct.ama1"))
print(rt.compare(s, "dct.ama1arf3"))

}


