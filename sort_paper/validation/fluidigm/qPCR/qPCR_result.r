# Analysis of the qPCR results.

source("git/utils.r")

data.dir = "data/expression/qPCR/20150428 daf-19 hlh-6/"

r1 = read.tsv(paste0(data.dir, "results1.tsv"))
r2 = read.tsv(paste0(data.dir, "results2.tsv"))

# utility to improve gene names
add.hyphen = function(s) {
  pattern = "^([a-z]+)([0-9]+)$"
  f = function(x) {
    a = regmatches(x, regexec(pattern, x))[[1]]
    if (length(a) > 0)
      x = paste0(a[2], "-", a[3])
    x
  }
  sapply(s, f)
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

# Utility to compute standard error.
se = function(x) {
  x = na.omit(x)
# XXX using sd() instead of sqrt(mean of squared diff)
#  s2 = mean( (x-mean(x)) ^ 2)
  sd(x) / sqrt(length(x))
}

# Compares one set of samples with N2 controls.
# Args:
#   s - a data.table of qPCR results
#   sample, target - the sample and target to get numbers for
# Returns: vector of ddct stats, and (unadjusted) p-value
rt.compare.1 = function(s, sample, target, dct.column) {
  dct.sample = s[ s$sample==sample & s$target==target , dct.column ]
  dct.N2 = s[ s$sample=="N2" & s$target==target , dct.column ]
  ddct = mean(dct.sample, na.rm=TRUE) - mean(dct.N2, na.rm=TRUE)
  ddct.se = se(dct.sample) + se(dct.N2)

  r = t.test(dct.sample, dct.N2, alternative="two.sided", var.equal=FALSE)
  data.frame(sample = sample, target = target,
    ddct = ddct, ddct.se = ddct.se,
    t = r$statistic, p = r$p.value,
    n.sample=length(dct.sample), n.N2=length(dct.N2), stringsAsFactors=FALSE)
}

# Compares all samples with N2 controls.
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

# Plots one set of statistics.
plot.stats = function(r, main) {
  ylim = range(c(2, r$ddct - r$ddct.se, r$ddct + r$ddct.se))
  par(mar=c(8,6,2,0) + 0.1)
  m = barplot(r$ddct, ylim=ylim, names.arg = r$target,
    col=c(rep("#ff0000a0", 3), rep("#0000ffa0", 3)),
    space=1, las=2, yaxt="n", main=main)
  mtext("          Relative expression", side=2, line=4.5) # XXX
  mtext("   daf-19 RNAi", side = 1, line=6, adj=0, col="#ff0000c0")
  mtext("hlh-6 RNAi   ", side = 1, line=6, adj=1, col="#0000ffc0")

  arrows(m, r$ddct - r$ddct.se, m, r$ddct + r$ddct.se,
    length=0.06, angle=90, code=3, col="black", lwd=2)
  y = trunc(ylim)[1] : trunc(ylim)[2]
  axis(2, at=y, labels = 2^y, las=1)

  text(m, 1.1, signif(r$p, 2), srt=90, cex=0.8)
}

# summarize the technical replicates
s1 = rt.summarize(r1)$small
# remove an outlier (even though it nominally worked)
s1 = s1[ !(s1$sample=="N2" & s1$rep=="022415c" & s1$target=="ama1") , ]
s2 = rt.summarize(r2)$small
# slight name adjustment
s2$rep = sub("022414c", "022415c", s2$rep)
s = rbind(s1, s2)

write.tsv(s, "git/sort_paper/validation/fluidigm/qPCR/qPCR_averages.tsv")

# compute dcts
s$ct.ama1 = matched.control(s, "ama1")
s$ct.arf3 = matched.control(s, "arf3")
s$dct.ama1 = s$ct.ama1 - s$ct
s$dct.ama1arf3 = (s$ct.ama1 + s$ct.arf3) / 2 - s$ct

# compute ddcts, and stats
rt.stats.1 = rt.compare(s, "dct.ama1")
print(rt.stats.1)
rt.stats = rt.compare(s, "dct.ama1arf3")
rt.stats$sample = add.hyphen(rt.stats$sample)
rt.stats$target = add.hyphen(rt.stats$target)
print(rt.stats)
write.tsv(rt.stats, "git/sort_paper/validation/fluidigm/qPCR/qPCR_result.tsv")

# subset this a bit
rt = rt.stats[ match(c("che-13", "C05D10.2", "ZK813.5",
  "phat-5", "F41G3.21", "T05B4.8"),
  rt.stats$target) , ]

# plot graph(s)
pdf("git/sort_paper/validation/fluidigm/qPCR/qPCR_result.pdf",
  width=3.5, height=6)
# plot.stats(rt.stats.1, "with ama-1 as control")
plot.stats(rt, "")
dev.off()

