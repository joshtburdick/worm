# Computes pairwise distance between two upstream motifs
# (as well as how often they occur in pairs.)

source("git/utils.r")

# running this on the filtered (non-redundant) motifs
motif.gene.dir = "git/cluster/motif/distAndConservation/5kb/"
known.motifs = {
  r = list.files(motif.gene.dir)
  sub("_upstreamMotifCons.tsv.gz", "", r)
}
motif.filter = read.tsv("git/tf/motif/motifFilter.tsv")
known.motifs.small =
  intersect(known.motifs, motif.filter$canonical.name)

# clustering to use
clustering1 = read.tsv(
  "git/cluster/hierarchical/hier.300.clusters/clusters.tsv")
clustering = clustering1[,2]
names(clustering) = rownames(clustering1)


# for gene names
x = read.tsv("git/cluster/readRatios.tsv")

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

# Upstream distances for all motifs.
upstream.dist = get.motif.upstream.dists(known.motifs.small)



# For each motif, gets the upstream distance of a random motif occurrence.
# rand.upstream.dist = function(m)
#   c(by(m$upstream.dist, m$gene, function(a) sample(a, 1)))


udt = upstream.dist.test(upstream.dist, clustering)
save(udt, file="git/sort_paper/tf/pairwise/pairwiseDist.Rdata")



