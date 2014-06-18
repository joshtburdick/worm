# Computes pairwise distance between two upstream motifs
# (as well as how often they occur in pairs.)

source("git/utils.r")
source("git/cluster/compute.cluster.centers.r")
source("git/data/name_convert.r")
source("git/sort_paper/tf/upstreamEnrichment.r")

# size of up-to-2kb-upstream-intergenic regions
upstream.size.2kb = apply(upstream.cons.dist[[2]], 1, sum)

# running this on the filtered (non-redundant) motifs
motif.gene.dir = "git/cluster/motif/distAndConservation/5kb/"
known.motifs = sub("_upstreamMotifCons.tsv.gz", "",
  list.files(motif.gene.dir))
motif.filter = read.tsv("git/tf/motif/motifFilter.tsv")
known.motifs.small =
  intersect(known.motifs, motif.filter$canonical.name)

# clustering to use
clustering1 = read.tsv(
  "git/cluster/hierarchical/hier.300.clusters/clusters.tsv")
clustering = clustering1[,2]
names(clustering) = rownames(clustering1)

# get cluster centers (and also gene names)
r = read.tsv("git/cluster/readRatios.tsv")
r.sort.only = r[ , 1:23 ]
r.center = compute.cluster.centers(r.sort.only, clustering)

# correlation of each cluster with known TFs
wtf = read.csv("data/tf/wTF2.1.csv", as.is=TRUE)
tf1 = unique(rename.gene.name.vector(
  union(wtf$Sequence.name.variant, wtf$Gene.public.name)))
# only consider genes which have at least some reads
tf2 = intersect(tf1, names(clustering))
tf.cluster.cor = cor(t(r.center), t(r.sort.only[tf2,]))

# Finds TFs correlated with a cluster center at some cutoff.
tf.cluster.cor.cutoff = function(tf.cluster.cor, r) {
  a = which(tf.cluster.cor >= r, arr.ind=TRUE)
  data.frame(cluster = as.character(rownames(tf.cluster.cor)[ a[,1] ]),
    tf = colnames(tf.cluster.cor)[ a[,2] ],
    stringsAsFactors = FALSE)
}

tfc = tf.cluster.cor.cutoff(tf.cluster.cor, 0.9)

# Gets the counts of one motif upstream of all genes.
get.motif.counts = function(m) {

  r = read.table(paste(motif.gene.dir, m, "_upstreamMotifCons.tsv.gz", sep=""),
    as.is=TRUE)
  colnames(r) =
    c("region.chr", "region.a", "region.b",
    "gene", "score", "strand",
    "motif", "motif.a", "motif.b", "motif.id", "motif.score",
    "motif.strand", "motif.cons")
  r$upstream.dist = ifelse(r$strand=="+",
    r$motif.b - r$region.b, r$region.a - r$motif.a)
  r
}

##### getting motifs per gene
# orthologs (using previous group of these rather than
# download from Wormbase, as it includes more genes, even
# though they may be not be so accurate)
motif.ortholog = read.table("git/tf/motif.ortholog.2.tsv",
  as.is=TRUE, header=TRUE)

# group similar motifs
motif.filter = read.tsv("git/tf/motif/motifFilter.tsv")

motif.ortholog$canonical.motif =
  motif.filter[ motif.ortholog$motif, "canonical.name" ]
motif.ortholog = motif.ortholog[ !is.na(motif.ortholog$canonical.motif) , ]

# picking one motif per gene, somewhat arbitrarily
motif.by.gene = {
  a = motif.ortholog
  a = a[ !duplicated(a$gene), ]
  b = motif.filter[ a$motif, "canonical.name" ]
  names(b) = a$gene
  b
}
# XXX
motif.by.gene = c(motif.by.gene, "pros-1" = "PROX1_DBD")
# 


most.significant.cutoffs =
  read.tsv("git/cluster/motif/enrichOptimize/mostSignificantCutoffs.hier.300.tsv")
rownames(most.significant.cutoffs) = most.significant.cutoffs$motif

# Gets upstream motif locations for the filtered set of motifs.
get.motif.upstream.dists = function(motif.list) {
  r = list()

  for(motif in motif.list) {
    write.status(motif)
    m = get.motif.counts(motif)
    m = m[ (m$upstream.dist >= -2000) &
      (m$motif.score >= most.significant.cutoffs[motif,"motif.score"]) , ]
    r[[motif]] = c(by(m$upstream.dist, m$gene, c))
  }

  r
}

# A random sample of distances for each pair of genes.
motif.dist.sample = function(m)
  sapply(m, function(a) sample(a, 1))

# Computes p-value of motifs co-occuring. This includes:
# - a chi-squared test of the motifs occuring more often
#   upstream of the genes in the cluster
# - a Wilcoxon of the distance between a random pair of
#   genes being less.
upstream.dist.test = function(upstream.dist, clustering) {
  r = NULL
  n.genes = length(clustering)

  for(m1 in names(upstream.dist))
    for(m2 in names(upstream.dist)) {
      r.1 = NULL
      if (m1 < m2) {

        # for each motif, pick one upstream location
        # ??? possibly move this inside the loop?
        dist.sample.1 = motif.dist.sample(upstream.dist[[m1]])
        dist.sample.2 = motif.dist.sample(upstream.dist[[m2]])
#print(dist.sample.1[1:10])
#print(dist.sample.2[1:10])
        for(cl in unique(sort(as.character(clustering)))) {
          write.status(paste(m1, m2, cl))

          cl.genenames = names(clustering)[ clustering==cl ]

          genes.with.both = intersect(names(dist.sample.1), names(dist.sample.2))
          genes.with.motif.1 = intersect(names(dist.sample.1), cl.genenames)
          genes.with.motif.2 = intersect(names(dist.sample.2), cl.genenames)

          # presumably if these occur in all the genes uniformly at random,
          # then this is the probability of both motifs being upstream of a gene
          p.both.indep = (length(dist.sample.1) / n.genes) *
            (length(dist.sample.2) / n.genes)

          names.both = intersect(genes.with.motif.1, genes.with.motif.2)
          num.in.cluster = length(names.both)

          # test for enrichment of co-occurrence near genes
# print(num.in.cluster)
# print(p.both.indep)
# print(cl.genenames)
# print(length(cl.genenames))
          enrich = num.in.cluster / (p.both.indep * length(cl.genenames))
          count.stat = NA
          count.p = 1

          if (enrich > 1) {
            tryCatch({
              ct = chisq.test(c(num.in.cluster, length(cl.genenames) - num.in.cluster),
                p = c(p.both.indep, 1 - p.both.indep), rescale.p=TRUE)
              count.stat = as.numeric(ct["statistic"])
              count.p = as.numeric(ct["p.value"])
            }, error = function(x) count.p = 1, warning = function(x) count.p = 1)
          }

          # test for motifs being more near by than expected
          dist.sample = abs(dist.sample.1[genes.with.both] - dist.sample.2[genes.with.both])
          dist.sample = dist.sample[ !is.na(dist.sample) ]

          dist.in.cluster = dist.sample[ intersect(names(dist.sample), names(clustering)[clustering==cl]) ]
          dist.not.in.cluster = dist.sample[ intersect(names(dist.sample), names(clustering)[clustering!=cl] ) ]
#          print(dist.in.cluster[1:3])
 #         print(dist.not.in.cluster[1:3])
# print(c(length(dist.in.cluster), length(dist.not.in.cluster)))

          dist.stat = NA
          dist.p = 1
#          tryCatch({
          if (length(dist.in.cluster) >= 1 & length(dist.not.in.cluster) >= 1) {
            w = wilcox.test(dist.in.cluster, dist.not.in.cluster, alternative="less")
            dist.stat = as.numeric(w$stat)
            dist.p = w$p.value
}
# },
#            error = function(x) 1, warning = function(x) 1)

          r.1 = rbind(r.1,
            data.frame(motif1 = m1, motif2 = m2, cluster = cl,
              n1 = length(genes.with.motif.1),
              n2 = length(genes.with.motif.2),
              n12 = length(names.both),
              enrich = enrich,
              count.stat = count.stat, count.p = count.p,
              mean.dist.cl = mean(dist.in.cluster),
              mean.dist.bg = mean(dist.not.in.cluster),
              dist.stat = dist.stat, dist.p=dist.p))
        }
      }
      # XXX hack to hopefully avoid slowness of rbind()
      r = rbind(r, r.1)
    }
  r
}

# Finds average difference between two vectors (as for
# upstream distance.)
# Args: dist.1, dist.2 - vectors of upstream distances
#   for two different genes
# Returns: average of their differences
avg.distance = function(dist.1, dist.2) {
  g = intersect(names(dist.1), names(dist.2))
  f = function(x) mean(abs(outer(dist.1[[x]], dist.2[[x]], "-")))
  sapply(g, f)
}

# Similar, but instead returns the "average closest"
# motif location (that is, for each motif A, the average
# distance to the closest motif B)
avg.closest.distance = function(dist.1, dist.2) {
  g = intersect(names(dist.1), names(dist.2))
  f = function(x) {
    a = abs(outer(dist.1[[x]], dist.2[[x]], "-"))
    mean(c(apply(a, 1, min), apply(a, 2, min)))
  }
  sapply(g, f)
}

# Returns the number of pairs of motifs which are within
# a small distance.
num.nearby.pairs = function(dist.1, dist.2) {
  g = intersect(names(dist.1), names(dist.2))
  f = function(x) {
    a = abs(outer(dist.1[[x]], dist.2[[x]], "-"))
    sum(a <= 200)
  }
  sapply(g, f)
}

# Tests for enrichment of pairs of motifs (possibly)
# corresponding to genes which are correlated with that cluster.
# Args:
#   clustering - a clustering (as a numeric vector, indexed by gene)
#   tf.cluster.cor - correlation of each tf with each cluster
#   cor.cutoff - the cutoff for correlation
#   upstream.dist - the upstream distances for all motifs
#   motif.by.gene - the motif to use for each gene
# Returns: data frame with columns:
#   cluster - the cluster
#   tf.1, tf.2 - the transcription factors
#   motif.1, motif.2 - the corresponding motifs
#   
correlated.tf.upstream.dist.test =
    function(clustering, tf.cluster.cor, cor.cutoff, upstream.dist, motif.by.gene) {
  r = NULL

  tfc = tf.cluster.cor.cutoff(tf.cluster.cor, cor.cutoff)

  for(cl in sort(unique(tfc$cluster))) {
    tfs = tfc[ tfc$cluster==cl, "tf" ]
    motifs = motif.by.gene[ tfs ]
    motifs = sort(unique(motifs[ !is.na(motifs) ]))

    for(m1 in motifs)
      for(m2 in motifs)
        if (m1<m2) {
          write.status(paste(cl, m1, m2))
          motif.dists = avg.closest.distance(
            upstream.dist[[m1]], upstream.dist[[m2]])

          # chi-squared test of the combined motif being enriched
          m = rep(0, length(clustering))
          names(m) = names(clustering)
          m[names(motif.dists)] = 1          
          motif.enrich = compute.enrichment.1(
            clustering == cl, m, upstream.size.2kb)

          # similarly, but for motifs being nearby
          m.nearby = rep(0, length(clustering))
          names(m) = names(clustering)
          nnp = num.nearby.pairs(upstream.dist[[m1]], upstream.dist[[m2]])     
          m.nearby[names(nnp)] = 1        
          motif.nearby.enrich = compute.enrichment.1(
            clustering == cl, m.nearby, upstream.size.2kb)

          # Wilcoxon of motifs being closer in the cluster
          dist.in.cluster = motif.dists[
            intersect(names(motif.dists), names(clustering)[clustering==cl]) ]
          dist.not.in.cluster = motif.dists[
            intersect(names(motif.dists), names(clustering)[clustering!=cl]) ]
          dist.stat = NA
          dist.p = 1
          if (length(dist.in.cluster) >= 1 & length(dist.not.in.cluster) >= 1) {
            w = wilcox.test(dist.in.cluster, dist.not.in.cluster)
            dist.stat = as.numeric(w$stat)
            dist.p = w$p.value
          }

          r = rbind(r, data.frame(
            cluster = cl,
            motif.1 = m1, motif.2 = m2,
            genes.1 = paste(intersect(names(motif.by.gene)[motif.by.gene==m1],
              tfs), collapse=","),
            genes.2 = paste(intersect(names(motif.by.gene)[motif.by.gene==m2],
              tfs), collapse=","),
            enrichment = motif.enrich["enrichment"],
            enrich.chisq = motif.enrich["chisq"],
            enrich.p = motif.enrich["p"],
            nearby.enrich = motif.nearby.enrich["enrichment"],
            nearby.chisq = motif.nearby.enrich["chisq"],
            nearby.p = motif.nearby.enrich["p"],
            dist.in.cluster = mean(dist.in.cluster, na.rm=TRUE),
            dist.not.in.cluster = mean(dist.not.in.cluster, na.rm=TRUE),
            dist.stat = dist.stat,
            dist.p = dist.p,
            stringsAsFactors=FALSE))
        }
  }

  rownames(r) = NULL
  r$enrich.p.corr = p.adjust(r$enrich.p, method="fdr")
  r$dist.p.corr = p.adjust(r$dist.p, method="fdr")
  r
}

# Upstream distances for all motifs.
# upstream.dist = get.motif.upstream.dists(known.motifs.small)
# save(upstream.dist, file="git/sort_paper/tf/pairwise/upstreamDist.Rdata")
# load("git/sort_paper/tf/pairwise/upstreamDist.Rdata")

# For each motif, gets the upstream distance of a random motif occurrence.
# rand.upstream.dist = function(m)
#   c(by(m$upstream.dist, m$gene, function(a) sample(a, 1)))


# udt = upstream.dist.test(upstream.dist, clustering)
# save(udt, file="git/sort_paper/tf/pairwise/pairwiseDist.Rdata")


for(cor.cutoff in c(0.9, 0.8, 0.7, 0.6, 0.5)) {
  cat(paste0(cor.cutoff, "\n"))
  r = correlated.tf.upstream.dist.test(clustering,
    tf.cluster.cor, cor.cutoff, upstream.dist, motif.by.gene) 
  write.tsv(r, paste0("git/sort_paper/tf/pairwise/pairwiseDist_", cor.cutoff, ".tsv"))
}

