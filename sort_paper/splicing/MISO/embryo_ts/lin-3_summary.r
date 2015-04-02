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
for(e in unique(lin3.events$event_name)) {
  r = lin3.events[ lin3.events$event_name == e , ]
  r = r[ order(r$t) , ]
  plot(r$t, r$miso_posterior_mean, type="b", pch=20,
    ylim=c(0,1), main=e)
  par(new=TRUE)
  plot(r$t, r$ci_low, type="b", pch=20,
    ylim=c(0,1), main=e, col="red")
  par(new=TRUE)
  plot(r$t, r$ci_high, type="b", pch=20,
    ylim=c(0,1), main=e, col="blue")
}
dev.off()

