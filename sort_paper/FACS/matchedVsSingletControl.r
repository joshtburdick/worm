# Enrichment of genes corresponding to markers used for sorting.

source("git/utils.r")
source("git/sort_paper/plot/experimentRename.r")

rpm = read.tsv("git/cluster/readsPerMillion.tsv")
# rpm.facs = rpm[ , c(1:21,36:54,65:68) ]

# pseudocount to add to read counts
pseudocount = 3

# log-transformed reads per million, with a pseudocount added
r = log2(pseudocount + rpm)
r = r[ , order(colnames(r)) ]

singlet.average =
  (r[,"cnd-1 singlets"] + r[,"pha-4 singlets"]) / 2

# Comparing read ratios using a matched control, versus the
# singlet control.
pos.genes = grep("^[^\\(]+\\(\\+\\)$", colnames(r), value=TRUE)
pos.genes = setdiff(pos.genes, c("hlh-16 (+)", "irx-1 (+)"))
neg.genes = sub("\\(\\+\\)", "(-)", pos.genes)

enrich = r[ , pos.genes ] - r[ , neg.genes ]
enrich.vs.singlet = r[ , pos.genes ] - singlet.average

png("git/sort_paper/FACS/matchedVsSingletControl.png",
  width=1400, height=1400)
par(mfrow=c(4,4))
par(mar=c(5,5,4,4)+0.1)
# XXX PDF file is huge
# pdf("git/sort_paper/FACS/matchedVsSingletControl.pdf",
#   width=7.5, height=10)
# par(mfrow=c(3,2))

fractions = setdiff(colnames(enrich), "pha-4 5/9 (+)")
fractions = fractions[ order(prettify.read.ratio.columns(fractions)) ]
for(j in fractions) {
  r = range(c(enrich[,j], enrich.vs.singlet[,j]))
  plot(enrich[,j], enrich.vs.singlet[,j],
    main=sub("\\(\\+\\) ", "", prettify.read.ratio.columns(j)),
    xlab="Enrichment", ylab="Enrichment vs. singlets",
    xlim=r, ylim=r, pch=20, cex=2, col="#00000040", lwd=2,
    cex.axis=1.5, cex.lab=2, cex.main=2)
  abline(a=0, b=1, col="#00000080", lwd=2)
}

dev.off()

