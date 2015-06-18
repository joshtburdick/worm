# Counts the number of things enriched in various
# clustered datasets.

source("git/utils.r")

# thresholds to try using
p.cutoffs = 10^c(-1:-15)

# Counts things enriched in one directory.
count.enriched = function(d) {
  r = NULL

  for(f in list.files(d)) {
    write.status(f)
    enrich = NULL
    load(paste0(d, "/", f))
    e1 = enrich
    p.most.sig = apply(e1[,,"p.corr",,,], c(1,2), min)

    # loop through significance cutoffs
    for(p.cutoff in p.cutoffs) {

      ps = p.most.sig <= p.cutoff
      write.status(paste(f, p.cutoff))

      # add to counts
      r1 = data.frame(name = sub(".Rdata$", "", f),
        p.cutoff = p.cutoff,
        num.zero = sum(p.most.sig == 0),
        total = sum(ps),
        num.motifs = sum(apply(ps, 1, any)),
        num.clusters = sum(apply(ps, 2, any)))
      r = rbind(r, r1)
    }
  }

  r
}

enriched.counts = count.enriched(
  "git/sort_paper/tf/motif/hyperg/hughes_20141202/")



