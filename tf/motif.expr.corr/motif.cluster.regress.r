# Looks for TFs, whose motifs predict membership in a cluster
# (where "membership" is defined as "correlation of expression
# with cluster mean")

# used to avoid cases in which motifs are too similar
source("git/tf/motif/motiv.utils.r")

load("git/tf/motif/meme.format.pwm.Rdata")

source("git/utils.r")

backspace.string = paste(rep("\b", 60), collapse="")

x = read.tsv("git/unmix/seq/cluster/readsFACSandTS.tsv")

# remove cases with no variance
x = x[ apply(x,1,var) > 0 , ]

# motif counts
load("git/tf/motif/motifCount/motif.counts.Rdata")
cat("loaded motifs\n")

# ChIP data
chip = as.matrix(read.tsv(gzfile("git/tf/chip/TF_chip_5kbUp.tsv.gz")))

# clusters
clusters = read.tsv("git/cluster/clusters.tsv")

# only keep genes which are in all lists
g = intersect(intersect(rownames(x), rownames(known.motifs)),
  rownames(clusters))
# ??? should this include all the genes?

x = x[g,]
motif = known.motifs[g,]
clusters = clusters[g,]
chip = chip[g,]

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

# Computes cluster centers, and correlations of each gene with
# each cluster center (which we consider "fuzzy cluster membership")
# Args:
#   expr - the expression data
#   clusters - the cluster centers of each gene
# Returns: a list with elements
#   cluster.center - centers of each cluster (currently mean)
#   cluster.corr - correlation of each gene with each cluster center
cluster.centers = function(expr, cluster) {

  # restrict to cases in which cluster is known
  expr = expr[ !is.na(cluster), ]
  cluster = cluster[ !is.na(cluster) ]

  # compute cluster centers
  cluster.center = matrix(NA, nrow=length(unique(cluster)), ncol=ncol(expr))
  rownames(cluster.center) = unique(cluster)
  colnames(cluster.center) = colnames(expr)

  for(cl in unique(cluster)) {
    cat(cl, "")
    cluster.center[cl,] =
      apply(expr[!is.na(cluster) & cluster==cl,], 2, mean)
  }

  list(cluster.center = cluster.center,
    cluster.corr = cor(t(expr), t(cluster.center)))
}

# Computes full linear model for one pair of TFs.
# Args:
#   ortho - table of orthologs
#   motif - matrix of motif counts upstream of each gene
#   cluster.corr - fuzzy cluster membership for each gene
#     (this should have named rows)
#   output.file - where to save results to
# Side effects: saves table in output.file (after each cluster)
# Returns: table of significant predictors of fuzzy membership
#   in each cluster (the same thing it writes to the file)
cluster.pairwise.interaction.regress =
    function(ortho, motif, cluster.corr, output.file) {

  r = NULL

  # for each cluster, and pair of motifs...
  for(cl in colnames(cluster.corr)) {
    for(i in 1:(nrow(ortho)-1)) {
      for(j in (i+1):nrow(ortho)) {

        # names of genes and motifs
        g1 = ortho[i,"gene"]
        g2 = ortho[j,"gene"]
        m1 = ortho[i,"motif"]
        m2 = ortho[j,"motif"]

        # find distance between motifs
        motif.dist = motifDistances(
          list(a=meme.format.pwm[[m1]], b=meme.format.pwm[[m2]]))[1]

        if (!(m1==m2) && motif.dist >= 0.05) {

          cat(backspace.string, cl, m1, m2, "  ")

          d = data.frame(y = cluster.corr[,cl], m1 = motif[,m1], m2 = motif[,m2])
          a = summary(lm(y ~ m1*m2, data=d))$coefficients
#        if("m1" %in% rownames(a)) {
          r = rbind(r, data.frame(cluster=cl, g1=g1, m1=m1, g2=g2, m2=m2,
            type="m1", t=a["m1",3], p=a["m1",4], p.adj=NA))
#        }
#        if("m2" %in% rownames(a)) {
          r = rbind(r, data.frame(cluster=cl, g1=g2, m1=m2, g2=g1, m2=m1,
            type="m1", t=a["m2",3], p=a["m2",4], p.adj=NA))
#        }
#        if("m1:m2" %in% rownames(a)) {
          r = rbind(r, data.frame(cluster=cl, g1=g1, m1=m1, g2=g2, m2=m2,
            type="m1:m2", t=a["m1:m2",3], p=a["m1:m2",4], p.adj=NA))
#          }
        }
      }
    }

#    r$p.adj = p.adjust(r$p, "fdr")

    # save file so far
    write.tsv(r, output.file)
  }

  r
}

# Regresses "membership on a cluster" (using the "fuzzy" definition
# of "correlation with cluster center") on some predictors.
# Args:
#   cluster.corr - correlations of genes with cluster centers
#   x - matrix of predictors
#   output.file - where to write output
# Side effects: writes to output file
cluster.corr.regress = function(cluster.corr, x, output.file) {
  cluster.corr.regress = NULL

  cat("\n")
  for(i in 1:ncol(cluster.corr)) {
    cat(backspace.string, i)
    y = cluster.corr[,i]
    ms = summary(lm(y ~ x))

    cluster.corr.regress[[i]] = sapply(c("call", "terms",
      "coefficients", "sigma", "df", "r.squared",
      "adj.r.squared", "fstatistic"), function(l) ms[[l]])
    if (log2(i) %% 1 == 0)
      save(cluster.corr.regress, file=output.file)
  }
  cat("\n")

  save(cluster.corr.regress, file=output.file)
}

# Regresses "membership in a cluster" on some predictors,
# using a GLM with a binomial link function.
# Args:
#   cl - a clustering
#   x - matrix of predictors
#   output.file - where to write output to
# Side effects: writes to output file (every so often, and
#   when all results have been computed)
cluster.membership.regress.glm = function(cl, x, output.file) {
  cluster.regress = NULL

  for(i in cl) {
    cat(backspace.string, i)
    y = cl == i
    ms = summary(glm(y ~ x, family = binomial(link = "logit")))
    m1 = sapply(c("call", "terms", "coefficients", "deviance", "aic",
      "df.residual", "null.deviance", "df.null", "iter",
      "dispersion", "df"), function(l) ms[[l]])

    cluster.regress[[i]] = m1
    if (log2(i) %% 1 == 0)
      save(cluster.regress, file=output.file)
  }

  save(cluster.regress, file=output.file)
}

cc = cluster.centers(x, clusters[,"hier.ts.100.clusters"])

# only keeping one gene with each motif
# XXX this results in only 178 motifs, which seems a bit low
ortho1 = ortho[ !duplicated(ortho$motif) , ]

# XXX for now, testing with just smattering of motifs
# ortho1 = ortho1[ sample(nrow(ortho1), 10), ]

# older model, including interaction terms
if (FALSE) {
cluster.pairwise.motif = cluster.pairwise.interaction.regress(
  ortho1, motif, cc$cluster.corr,
  "git/tf/motif.expr.corr/motif.cluster.regress.ts.100.tsv")
# write.tsv(cluster.pairwise.motif, "git/tf/motif.expr.corr/motif.cluster.regress.ts.100.tsv")
}

# GLM model
if (FALSE) {
system("mkdir -p git/tf/motif.expr.corr/glm/")
cluster.membership.regress.glm(clusters[,"hier.ts.100.clusters"], motif,
  "git/tf/motif.expr.corr/glm/hier.ts.100.clusters_motif.Rdata")
}

# linear model
if (TRUE) {
  system("mkdir -p git/tf/motif.expr.corr/cluster.lm/")
  cluster.corr.regress(cc$cluster.corr, chip,
    "git/tf/motif.expr.corr/cluster.lm/hier.ts.100.clusters_chip.Rdata")
  cluster.corr.regress(cc$cluster.corr, cbind(motif, chip),
    "git/tf/motif.expr.corr/cluster.lm/hier.ts.100.clusters_motif_chip.Rdata")
  cluster.corr.regress(cc$cluster.corr, motif,
    "git/tf/motif.expr.corr/cluster.lm/hier.ts.100.clusters_motif.Rdata")

}

