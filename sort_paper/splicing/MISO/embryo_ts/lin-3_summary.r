# Summarizes just the events relevant to lin-3.

source("git/utils.r")

miso.all = read.tsv(gzfile("git/sort_paper/splicing/MISO/embryo_ts/misoSummary.tsv.gz"))

# process annotation somewhat
miso.all$t = as.numeric(sub(".*Early-Embryos-50-", "", miso.all$file))

miso.all$mrna_start = sapply(strsplit(miso.all$"mRNA_starts", ","), function(a) a[1])

# select just things in the lin-3 region, near
# IV:11,060,174-11,061,494
lin3.events = miso.all[ miso.all$chrom=="IV" &
  miso.all$mrna_start >= 11060174 &
  miso.all$mrna_start <= 11061494 , ]

pdf("git/sort_paper/splicing/MISO/embryo_ts/lin-3_summary.pdf",
  width=7.5, height=10)
par(mfrow=c(4,1))
par(mar=c(4,4,4,1)+0.1)
a = unique(lin3.events[ , c("as.type", "event_name") ])

for(i in 1:nrow(a)) {
  e = a[ i , "event_name" ]
  r = lin3.events[ lin3.events$event_name == e , ]
  r = r[ order(r$t) , ]

  plot(r$t, r$miso_posterior_mean, type="l", pch=20,
    ylim=c(0,1),
    xlab="Time (minutes)", ylab="Proportion spliced in",
    main = a[ i, "as.type" ])
  mtext(a[ i, "event_name"], line = 0.5, cex=0.7)
  arrows(r$t, r$ci_low, r$t, r$ci_high, angle=90,
    length=0.05, code=3, col="#00000080")
}
dev.off()

