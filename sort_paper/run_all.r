# R script to re-run, if not all, then at least a
# substantial part of the analysis.
# FIXME: automate dependency checking?

source("git/unmix/seq/timing/deconvolved_embryo_ts.r")
source("git/cluster/readRatios.r")

source("git/cluster/hierarchical.r")

source("git/data/wormbase/anatomy.ontology.group.r")
source("git/data/wormbase/expr.cluster.r")

# freezing this for now
# source("git/data/homology/wormbase.tf.ortholog.r")

source("git/tf/motif/meme.tf.annotate.r")

source("git/sort_paper/enrichment/anatomyEnrichment.r")


