# Compares FACS-sorted fractions with whole-embryo data.
# This is partly as a way of QCing, and also is related
# to unmixing.

source("git/utils.r")

rpm = read.tsv("git/cluster/readsPerMillion.tsv")
rpm$ungated.mean = (rpm$"pha-4 ungated" + rpm$"cnd-1 ungated") / 2
rpm$singlet.mean = (rpm$"pha-4 singlet" + rpm$"cnd-1 singlet") / 2

# Regresses whole-embryo data on flow-sorted fractions,
# and plots the result.
regress.on.facs = function(sort.fraction) {
  f = function(x) { x[x<0] = 0; log2(1+x) }
  r1 = data.frame(ungated = rpm[ , "ungated.mean" ],
    singlet = rpm[ , "singlet.mean" ],
    plus = rpm[ , paste0(sort.fraction, " (+)") ],
    minus = rpm[ , paste0(sort.fraction, " (-)") ])

  m1 = lm(ungated ~ plus + minus + 0, data = r1)
  m2 = lm(ungated ~ plus + minus + singlet + 0, data = r1)

  plot(f(r1$ungated), f(predict(m1)),
#  plot(log2(1+r1$ungated), log2(1+predict(m1)),
    main = sort.fraction,
    xlab="whole embryo rpm", ylab="regressed rpm",
    pch=20, col="#00000080")
  plot(f(r1$ungated), f(predict(m2)),
#  plot(log2(1+r1$ungated), log2(1+predict(m2)),
    main = paste(sort.fraction, "including singlets"),
    xlab="whole embryo rpm", ylab="regressed rpm",
    pch=20, col="#00000080")
  c(cor.orig = cor((r1$ungated), (predict(m1))),
    cor.w.singlets = cor((r1$ungated), (predict(m2))))
}



png.dir = "git/sort_paper/FACS/regress/regressionPlots"

system(paste0("mkdir -p ", png.dir))
r = NULL
for(sf in c("ceh-27", "ceh-36", "ceh-6",
    "cnd-1 12/14", "cnd-1 1/4", "cnd-1 8/19",
    "F21D5.9", "mir-57", "mls-2", "pal-1",
    "pha-4 12/9", "pha-4 5/9", "pha-4 9/1",
    "pros-1", "ttx-3", "unc-130")) {
  sf.name = gsub("/", "_", sf)
  png(paste0(png.dir, "/", sf.name, ".png"), width=800, height=400)
  par(mfrow=c(1,2)) 
  r1 = regress.on.facs(sf)
  r = rbind(r, r1)
  dev.off()
}

pdf("git/sort_paper/FACS/regress/facsRegress.pdf")
plot(r[,1], r[,2],
  pch=20, col="#00000080", lwd=2, xlim=c(0.5,1), ylim=c(0.5, 1),
  xlab="Correlation, using (+) and (-)",
  ylab="Correlation, using (+), (-), and singlets")
abline(0,1, col="#00000080", lwd=2)
dev.off()

print(wilcox.test(r[,1], r[,2], alternative="less", paired=TRUE))


