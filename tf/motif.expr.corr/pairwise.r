# Looks for pairs of correlated TFs, in which presence of
# upstream motifs significantly predicts correlation with one of those TFs.

# used to avoid cases in which motifs are too similar
source("git/tf/motif/motiv.utils.r")

load("git/tf/motif/meme.format.pwm.Rdata")

source("git/utils.r")

# x = read.tsv("git/unmix/seq/cluster/readsFACSandTS.tsv")
x = read.tsv("git/cluster/readRatios.tsv")

# XXX: omitting timeseries data for now
# x = x[,1:23]

# Computes a weighted correlation, based on
# http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient#Calculating_a_weighted_correlation
# Args:
#   x, y - matrices to compute correlations of
#   w - non-negative weights
#     (with length equal to the number of rows in x and y)
# Returns: weighted correlations of columns of x and y
weighted.cor = function(x, y, w) {
  gc()

  # weighted mean of each column of x
  w.mean = function(x, w)
    apply(w * x, 2, sum) / sum(w)

  # weighted covariance
  w.cov = function(x, y, w)
    t(w * t(t(x) - w.mean(x,w))) %*% t(t(y) - w.mean(y,w)) / sum(w)

  w.cov(x, y, w) / sqrt( diag(w.cov(x, x, w)) %o% diag(w.cov(y, y, w)) )
}

# remove cases with no variance
x = x[ apply(x,1,var) > 0 , ]

# motif counts
load("git/tf/motif/motifCount/motif.counts.Rdata")
cat("loaded motifs\n")

# only keep genes which are in both lists
g = intersect(rownames(x), rownames(known.motifs))
x = x[g,]
motif = known.motifs[g,]

# instances of orthologs
ortho = read.table("git/tf/motif.ortholog.2.tsv", header=TRUE, as.is=TRUE)

# XXX check for other gene names?
ortho = ortho[ ortho$gene %in% rownames(x), ]

# just keep cases in which there's a motif, and only keep
# the motif from the "best" homolog (by "gene.score"; arguably
# these should be normalized in some way)
ortho = ortho[ ortho$motif != "" , ]
ortho = ortho[
  order(ortho$gene,
    ortho$species,
    -ortho$related.gene.score) , ]

# XXX some of these motifs aren't included (I'm not sure why),
# so omit those cases
ortho = ortho[ ortho$motif %in% colnames(motif) , ]

# ortho = ortho[ !duplicated(ortho[ , c("gene", "species"),]) , ]
# ??? just keeping unique gene names
ortho = ortho[ !duplicated(ortho[ , c("gene"),]) , ]

cat("filtered orthologs\n")

tfs = intersect(rownames(x), ortho$gene)

# how much each gene is correlated with those transcription factors
tf.cor = cor(t(x), t(x[tfs,]))
# tf.cor = weighted.cor(t(x), t(x[tfs,]), c(rep(1,23), rep(0.5,15)))

cat("computed correlations\n")

# Computes full linear model for one pair of TFs.
# Args:
#   ortho - data frame including columns "gene" and "motif", indicating
#     instances in which a gene 
#   motif - matrix of motif counts upstream of each gene
#   tf.cor - matrix of correlation of each gene with each TF
#   min.cor - correlation between TF_i and TF_j must be at least this high
# Returns: regression object corresponding to multiple regression of
# correlation of genes with TF_i, regressed on number of motifs of
# TF_i and TF_j.
pairwise.interaction.regress = function(ortho, motif, tf.cor, min.cor) {
  gc()

  r = NULL

  for(i in 1:(nrow(ortho)-1)) {
    for(j in (i+1):nrow(ortho)) {
#  for(i in 19) {
#    for(j in 176) {
      g1 = ortho[i, "gene"]
      g2 = ortho[j, "gene"]

      m1 = ortho[i, "motif"]
      m2 = ortho[j, "motif"]

      # distance between these motifs
      motif.dist = motifDistances(
        list(a=meme.format.pwm[[m1]], b=meme.format.pwm[[m2]]))[1]

      if (!(m1==m2) && (tf.cor[g1,g2] >= min.cor) &&
        motif.dist >= 0.05) {
cat(i, j, "  ")

        d = data.frame(y = tf.cor[,i], m1 = motif[,m1], m2 = motif[,m2])
        a = summary(lm(y ~ m1*m2, data=d))$coefficients
#        if("m1" %in% rownames(a)) {
          r = rbind(r, data.frame(g1=g1, g2=g2, type="m1", t=a["m1",3], p=a["m1",4]))
#        }
#        if("m2" %in% rownames(a)) {
          r = rbind(r, data.frame(g1=g2, g2=g1, type="m1", t=a["m2",3], p=a["m2",4]))
#        }
#        if("m1:m2" %in% rownames(a)) {
          r = rbind(r, data.frame(g1=g1, g2=g2, type="m1:m2", t=a["m1:m2",3], p=a["m1:m2",4]))
#        }
      }
    }
  }

  r$p.adj = p.adjust(r$p)
  r
}


# Utility to convert two vectors of strings to a string pair,
# separated by some character.
sorted.string.pair = function(a, b, sep="_") {
  a1 = as.character(a)
  b1 = as.character(b)
  paste(pmin(a1,b1), pmax(a1,b1), sep=sep)
}

# Counts how many interactions were "filtered" by this analysis,
# compared to just using this as a coexpression network.
count.pairwise.filtered = function() {

  # first, count TFs correlated at least this much
  tf1 = tf.cor[colnames(tf.cor),]
  cat("num pairs of TFs correlated >= 0.7 =",
    (sum(tf1 >= 0.7) - 462) / 2, "\n")

  # then, count the number included by this regression analysis
  load("git/tf/motif.expr.corr/pairwise.Rdata")
  r1 = pairwise.0.7
  r1 = r1[ r1$p.adj <= 0.05 , ]
  cat("number with any significant regression =",
    length(unique(sorted.string.pair(r1$g1, r1$g2))), "\n")

  # number of interaction terms
  cat("number of interaction terms = ", sum(r1$type=="m1:m2"), "\n")



}

if (FALSE) {
  pairwise.0.9 = pairwise.interaction.regress(ortho, motif, tf.cor, 0.9)
  pairwise.0.8 = pairwise.interaction.regress(ortho, motif, tf.cor, 0.8)
  pairwise.0.7 = pairwise.interaction.regress(ortho, motif, tf.cor, 0.7)

  save(pairwise.0.9, pairwise.0.8, pairwise.0.7,
    file="git/tf/motif.expr.corr/pairwise.Rdata")
}


