# R script to re-run, if not all, then at least
# some of the analysis.

# function to run a script (used to ensure that
# we only use the side effects of the script.)
run = function(script) {
  system(paste0("R --no-save CMD BATCH ", script))
}

run("git/unmix/seq/timing/deconvolved_embryo_ts.r")
run("git/cluster/readRatios.r")

run("git/cluster/hierarchical.r")



source("git/data/wormbase/anatomy.ontology.group.r")
source("git/data/wormbase/expr.cluster.r")

# freezing this for now
# source("git/data/homology/wormbase.tf.ortholog.r")

source("git/tf/motif/meme.tf.annotate.r")

source("git/sort_paper/enrichment/anatomyEnrichment.r")


source("git/sort_paper/enrichment/summary.r")


source("git/sort_paper/enrichment/clustersWithAnatomyAnnotation.r")


source("git/sort_paper/plot/web/clusters.r")

source("git/sort_paper/network/clusterTF.r")

source("git/sort_paper/plot/enrichment/stackedPlots.r")

source("git/cluster/motif/plot/distConservationWilcoxon.r")

# comparison of "number of motifs enriched for different sizes of cluster"
source("git/sort_paper/tf/motif/hyperg/numEnriched/clusteringComparison.r")
source("git/sort_paper/tf/motif/hyperg/numEnriched/numEnrichedInClustering.r")




