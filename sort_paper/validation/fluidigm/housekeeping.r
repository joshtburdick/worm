# Finds housekeeping genes.

source("git/utils.r")

# one of the clusterings
clustering = {
  cl1 = read.table("git/cluster/hierarchical/hier.300.clusters/clusters.tsv",
    sep="\t", header=TRUE, as.is=TRUE)
  cl = cl1[,3]
  names(cl) = cl1[,2]
  cl
}

# list of genes: those used in this clustering, which
# (presumably) are at least somewhat annotated
genes = grep("-", names(clustering), value=TRUE)

# reads per million
rpm = read.tsv("git/cluster/readsPerMillion.tsv")
rpm = rpm[ genes, !grepl("^(HS|RNAi)", colnames(rpm)) ]
lrpm = log10(1 + rpm)

# the enrichments
r = as.matrix(read.tsv("git/cluster/readRatios.tsv"))
r = r[ genes , ]
r = r[ , c(1:23) ]

# summary statistics for each gene
mean.lrpm = apply(lrpm, 1, mean)
sd.enrich = apply(r, 1, sd)

# restrict to just genes with at least some mean RPM
genes = names(mean.lrpm)[ mean.lrpm >= 2 ]
mean.lrpm = mean.lrpm[genes]
sd.enrich = sd.enrich[genes]

# simple ranking based on these (smaller is better)
# h.score = order(-mean.rpm) + 1000 * order(sd.enrich)
h.score = sd.enrich # / log10(1+mean.rpm)

names(h.score) = names(mean.lrpm)
h.score = sort(h.score)
n = 30
g1 = names(h.score[1:n])

# XXX
g1 = c(g1, "ama-1")

g1 = g1[ order(mean.lrpm[g1]) ]

# statistics for just those genes
r1 = r[ g1, ]
lrpm1 = lrpm[ g1, ]
mean.lrpm1 = apply(lrpm1, 1, mean)

h = rep(c(0,1/3,2/3), 100)[1:n]

pdf("git/sort_paper/validation/fluidigm/housekeeping.pdf",
  width=15, height=6)
xlim = range(mean.lrpm1)
enrich.max = max(abs(r1)) + 0.5
ylim = c(-enrich.max, enrich.max)
plot(1, 1, type="n", xlim=xlim, ylim=ylim,
  main="", xlab="Mean log10(1 + RPM)", ylab="Enrichment")
for(i in 1:n) {
  par(new=TRUE)
  plot(rep(mean.lrpm1[i], ncol(r1)), r1[ i, ],
    xlim=xlim, ylim=ylim,
    xlab="", ylab="", xaxt="n", yaxt="n",
    pch=20, col=hsv(h[i], 1, 0.8, 0.6))
}

abline(h=0, col="#00000060", lwd=2)
text(mean.lrpm1, apply(r1, 1, max),
  paste0("  ", g1), adj=c(0,0.5), col=hsv(h, 1, 0.8),
  cex=0.8, srt=90)
dev.off()

