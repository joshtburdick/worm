# Looks for anatomy terms which are overrepresented in one fraction or another.

wb.anatomy = read.table(
  gzfile("data/wormbase/anatomyTerm.ws220.tsv.gz"),
  header=TRUE, sep="\t", quote="", as.is=TRUE)

anatomy.term.to.name = by(wb.anatomy$Anatomy.Term, wb.anatomy$Anatomy.Term.ID,
  function(x) as.vector(x)[1])

r = as.matrix(read.table(gzfile(
    "git/unmix/seq/quant/readsPerMillion_pooled.tsv.gz"),
  header=TRUE, row.names=1, check.names=FALSE, as.is=TRUE))

lr = log2(1 + r)
le = lr - lr[,"all"]

if (FALSE) {
  pha4.plus.enriched.anatomy.terms = test.enriched(pha4.plus.enriched.genes)
  pha4.minus.enriched.anatomy.terms = test.enriched(pha4.minus.enriched.genes)
}

# Tests for enriched tissues, using a paired t-test.
enriched.t.test = function(expr.1, expr.2) {

  # restrict anatomy table to just genes we know
  genes = intersect(names(expr.1), names(expr.2))
#  print(genes[1:10])
  a = wb.anatomy[ wb.anatomy$Gene.Public.Name %in% genes , ]
#  print(dim(a))

  # Tests one anatomy term.
  t.test.1 = function(at) {
    g = a[ a$Anatomy.Term.ID == at , "Gene.Public.Name" ]
 
    x1 = expr.1[g]
    x2 = expr.2[g]

    # make sure enough data is present, and has some variance
    if (length(g) < 2 || (var(expr.1[g]) + var(expr.2[g]) <= 1e-3))
      return(NULL)

    # XXX unsure about "equal variances" assumption here
#    test = t.test(expr.1[g], expr.2[g], paired=TRUE, var.equal=TRUE)
    test = wilcox.test(expr.1[g], expr.2[g], paired=TRUE)
    data.frame(name = anatomy.term.to.name[at],
      statistic = test$statistic, n = length(g),
      p = test$p.value)    # formerly had "parameter = test$parameter"
  }

  r = NULL
  for(at in unique(wb.anatomy$Anatomy.Term.ID)) {
#    cat(at, "")
    test = t.test.1(at)
    if (!is.null(test)) {
      r = rbind(r, t.test.1(at))
    }
  }
#  cat("\n")
  r$p.corr = p.adjust(r$p, method="BH")

  r[ order(r$statistic), ]
}

# computes anatomy terms enriched in each fraction
compute.anatomy.enriched = function() {
  system("mkdir -p git/unmix/seq/cluster/anatomy_enriched_wilcox/")
  for(g in setdiff(colnames(lr), "all")) {
    cat(g, "")
    r = enriched.t.test(lr[,g], lr[,"all"])
    write.table(r[ r$p.corr <= 0.05, ], sep="\t", col.names=NA,
      file=paste("git/unmix/seq/cluster/anatomy_enriched_wilcox/", g, ".tsv", sep=""))
  }
}


