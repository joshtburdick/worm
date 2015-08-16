# Summarizes results of searches for enriched
# anatomy terms, and other Wormbase cluster terms.

source("git/utils.r")
source("git/sort_paper/enrichment/asTable.r")

# number of genes assumed expressed (used for computing fold-enrichment)
# this number is from "git/sort_paper/plot/numEnrichedInFractions.r"
num.genes.expressed = 15683

# various names of things
# (arguably these should be stored elsewhere)
load("git/data/wormbase/anatomy.ontology.group.Rdata")
anatomy.term.to.name =
  by(ao.group$group.name, ao.group$group, function(a) as.character(a)[1])

load("git/data/wormbase/expr.cluster.Rdata")
expr.cluster.name = expr.cluster.descr

# function for annotating phenotypes (RNAi and otherwise)
pheno.name = {
  r = read.tsv(gzfile("data/wormbase/phenotype_WS220.tsv.gz"))
  pheno.name = r[,1]
  names(pheno.name) = rownames(r)
  pheno.name
}
phenotype.annotate = function(p) {
  p1 = sub("RNAi_", "", p)
  ph = pheno.name[ p1 ]
  i = grepl("RNAi_", p)
  ph[ i ] = paste(ph[i], "(RNAi)")
  ph
}

# Summarizes just the significant results in a directory.
summarize.enrich = function(a) {
  r = NULL

  for(i in 1:dim(a)[1])
    for(j in 1:dim(a)[2])
      if (a[i,j,"p.corr"] <= 0.05) {
        a1 = a[i,j,]
        r = rbind(r,
          data.frame(group = dimnames(a)[[1]][i],
            cluster = dimnames(a)[[2]][j],
            num.in.group = a1["num.in.group"],
            num.in.cluster = a1["num.in.cluster"],
            num.intersect = a1["num.intersect"],
            enrichment = round((a1["num.intersect"] / a1["num.in.cluster"]) /
              (a1["num.in.group"] / num.genes.expressed), 2),
            p = signif(a1["p"], 3),
            p.corr = signif(a1["p.corr"], 3),
            stringsAsFactors=FALSE))
      }
  r
}

summarize.enrich.dir = function(in.dir, out.dir, annotate.f) {
  system(paste("mkdir -p", out.dir))

  for(f in list.files(in.dir)) {
    # was "hier.300.clusters"
    write.status(f)

    output.name = sub(".Rdata$", "", f);

    r = NULL
    load(paste0(in.dir, "/", f))
    s = summarize.enrich(r)

    s$group.name = annotate.f(s$group)

    s = s[ , c("cluster", "group", "group.name", "num.in.group",
      "num.in.cluster", "num.intersect", "enrichment", "p", "p.corr") ]
    s = s[ order(s$cluster, s$p.corr) , ]

    rownames(s) = NULL
    write.tsv(s, paste0(out.dir, "/", output.name, ".tsv"))
  }
}

# Writes out a summary containing all the data
# (not filtered by significance.)
summarize.all = function() {
  rd = "git/sort_paper/enrichment/"
  out.dir = "git/sort_paper/enrichment/summary_all/"

  load(paste0(rd, "anatomyEnrichment/facs.Rdata"))
  a = enrich.result.to.table(r)
  a$name = anatomy.term.to.name[ a$id ]
  a = a[,c(8,1:7)]
  write.tsv(a, paste0(out.dir, "anatomyEnrichment.facs.tsv"))

  load(paste0(rd, "anatomyEnrichment/hier.300.clusters.Rdata"))
  a = enrich.result.to.table(r)
  a$name = anatomy.term.to.name[ a$id ]
  a = a[,c(8,1:7)]
  write.tsv(a, paste0(out.dir, "anatomyEnrichment.hier.300.clusters.tsv"))

  load(paste0(rd, "wormbaseCluster/facs.Rdata"))
  a = enrich.result.to.table(r)
  a$name = expr.cluster.name[ a$id ]
  a = a[,c(8,1:7)]
  write.tsv(a, paste0(out.dir, "wormbaseCluster.facs.tsv"))

  load(paste0(rd, "wormbaseCluster/hier.300.clusters.Rdata"))
  a = enrich.result.to.table(r)
  a$name = expr.cluster.name[ a$id ]
  a = a[,c(8,1:7)]
  write.tsv(a, paste0(out.dir, "wormbaseCluster.hier.300.clusters.tsv"))
}

if (TRUE) {
summarize.enrich.dir("git/sort_paper/enrichment_cutoff1/anatomyEnrichment/",
  "git/sort_paper/enrichment/summary_cutoff1/anatomyEnrichment/",
  function(a) anatomy.term.to.name[a])

summarize.enrich.dir("git/sort_paper/enrichment_cutoff1/wormbaseCluster/",
  "git/sort_paper/enrichment/summary_cutoff1/wormbaseCluster/",
  function(a) expr.cluster.name[a])
}

# summarize.enrich.dir("git/sort_paper/enrichment/phenotype/allResults/",
#   "git/sort_paper/enrichment/summary/phenotype/",
#   phenotype.annotate)

# summarize.all()

