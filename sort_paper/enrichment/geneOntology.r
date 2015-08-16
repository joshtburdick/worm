# Annotate clusters with Gene Ontology info.

source("git/utils.r")
source("git/data/gostats.r")
source("git/sort_paper/enrichment/groupToList.r")
source("git/sort_paper/FACS/enrichedInFraction.r")

# for background genes
source("git/sort_paper/plot/numEnrichedInFractions.r")

output.dir = "git/sort_paper/enrichment/geneOntology/"
system(paste("mkdir -p", output.dir))

# determine "background" gene symbols
background.genes = names(max.expr)[ max.expr > 1 ]

# Annotates many sets of genes.
# Args:
#   background.genes - the background genes
#   gene - a vector of genes
#   cluster - the names of a group of genes
#     (note that a gene can be in more than one cluster)
# Returns: data frame with
#   set - name of the set of genes
#   ... enriched Gene Ontology info
go.annotate = function(background.genes, gene, cluster) {
  r = NULL

  for(cl in sort(unique(cluster))) {
    gh = go.hyperg(background.genes, gene[ cluster == cl ], c(0, 50))
    if (!is.null(gh)) {
      r1 = cbind(cluster=cl, gh)
#      print(dim(r1))
      r = rbind(r, r1)
    }
  }

  if (!is.null(r)) {
    r$p.corr = p.adjust(r$Pvalue, method="fdr")
    r = r[ r$p.corr <= 0.05 , ]
  }

  r
}

# Annotates many sets of genes (as above, except that this expects
# genes in a different format.)
# Args:
#   clusters - a list of boolean vectors
# Returns: a data.frame of results
go.annotate.list = function(clusters) {
  r = NULL
  num.tests = 0
  for(cl in names(clusters)) {
    write.status(cl)
    a = clusters[[cl]]
    if (sum(a) > 0) {
    gh = go.hyperg(names(a), names(a)[a], c(0, 50))
      if (!is.null(gh) && !is.null(gh$r) && nrow(gh$r) > 0) {
        r1 = cbind(cluster=cl, gh$r)
        r = rbind(r, r1)
        num.tests = num.tests + gh$num.tests
      }
    }
  }

  if (!is.null(r)) {
    r$p.corr = p.adjust(r$Pvalue, method="fdr", n=num.tests)
    r = r[ r$p.corr <= 0.05 , ]
  }

  r
}

lineage.enriched.go = function() {
  le = read.tsv("git/sort_paper/unmix/lineageEnriched.tsv")
  # le = le[1:60,]  # XXX testing
  e = go.annotate(background.genes, le$gene, le$set)
  write.tsv(e, paste0(output.dir, "lineageEnriched.tsv"))
}

cluster.enriched.go = function() {
  clustering.dir = "git/cluster/hierarchical/"
  for(f in c("hier.300.clusters")) {   # , list.files(clustering.dir))) {
    cat(f, "\n")
    clusters = cluster.to.gene.list(paste0(clustering.dir, f, "/clusters.tsv"))
    r = go.annotate.list(clusters)
    write.tsv(r, paste0(output.dir, "/", f, ".tsv"))
  }
}

facs.enriched.go = function() {
  # XXX for testing
#  a = list("singlets enriched" = facs.enriched.depleted[["singlets enriched"]],
#    "pha-4 enriched" = facs.enriched.depleted[["pha-4 enriched"]])
#  r = go.annotate.list(a)

  r = go.annotate.list(facs.enriched.depleted)
  write.tsv(r, paste0(output.dir, "/facs.tsv"))
}

# Checks for GO annotation in FACS-sorted samples at a looser
# threshold.
# Args:
#   cutoff - an enrichment cutoff
# Side effects: writes output (with the cutoff in the name)
facs.enriched.go.looser = function(cutoff) {
  r = go.annotate.list(facs.enriched.depleted)
  write.tsv(r, paste0(output.dir, "/facs_cutoff=", cutoff, ".tsv"))

}

# facs.enriched.go()
# cluster.enriched.go()


