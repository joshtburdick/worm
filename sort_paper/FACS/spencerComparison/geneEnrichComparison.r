# Compares genes enriched in sort fractions with
# genes enriched in Spencer data.

source("git/utils.r")
source("git/data/name_convert.r")
source("git/plot/heatmap.r")

# which cells are in which sort fraction / Spencer experiment
source("git/sort_paper/unmix/sortMatrix.r")

spencer.m = as.matrix(read.tsv(
  "git/sort_paper/validation/spencerSortMatrix.tsv"))
spencer.m = spencer.m[ sort(rownames(spencer.m)) , ]

source("git/sort_paper/FACS/enrichedInFraction.r")
# using a lower cutoff for "enriched" or "depleted"
ed1 = get.enriched.and.depleted(r.sort.only.averaged, log2(4))

# convert to lists of genes (as opposed to boolean vectors)
facs.1 = sapply(ed1, function(a) names(a)[a])
facs.1 = facs.1[ grep("enriched", names(facs.1), value=TRUE) ]
facs.1 = facs.1[ grep("singlets", names(facs.1), invert=TRUE, value=TRUE) ]
facs.1 = facs.1[ grep("ceh-6 \\(-\\) hlh-16 \\(-\\)",
  names(facs.1), invert=TRUE, value=TRUE) ]
names(facs.1) = sub(" enriched", "", names(facs.1))
facs.1 = facs.1[ sort(names(facs.1)) ]

m.unnormalized = m.unnormalized[ names(facs.1) , ]

# information about expression clusters from WormBase
load("git/data/wormbase/expr.cluster.Rdata")
wb.cluster = expr.cluster
wb.cluster$gene = rename.gene.name.vector(wb.cluster$gene)
wb.cluster = unique(wb.cluster)
wb.cluster$group.name = expr.cluster.descr[ wb.cluster$group ]

# just limiting to Spencer "selectively enriched" clusters for now
wb.cluster = wb.cluster[ , c("gene", "group") ]
wb.cluster = wb.cluster[
  grep("WBPaper00037950.*_embryo_SelectivelyEnriched", wb.cluster$group) , ]

# munge names
wb.cluster$group = sub("WBPaper00037950:", "", wb.cluster$group)
wb.cluster$group = sub("_embryo_SelectivelyEnriched", "", wb.cluster$group)

spencer.cluster = list()
for(g in sort(unique(wb.cluster$group)))
  spencer.cluster[[g]] = unique(wb.cluster[wb.cluster$group==g, "gene"])
# XXX match up names
spencer.cluster = spencer.cluster[c(1,3:12,2,13)]
names(spencer.cluster) = rownames(spencer.m)

# Converts a boolean matrix to a list of TRUE elements.
bool.matrix.to.list = function(a) {
  r = NULL
  for(i in rownames(a))
    r[[ i ]] = colnames(a)[ a[i,] ]
  r
}

# Jaccard similarity.
# Args: a, b - vectors
# Returns: Jaccard similarity of a and b.
jaccard.sim = function(a, b)
  length(intersect(a,b)) / length(union(a,b))

# Makes a matrix of how much things overlap.
# Args: a, b - two lists of sets of things
# Returns: matrix of proportions overlapping
overlap.matrix = function(a, b) {
  r = matrix(nrow=length(names(a)), ncol=length(names(b)))
  rownames(r) = names(a)
  colnames(r) = names(b)
  for(i in rownames(r))
    for(j in colnames(r))
      r[i,j] = jaccard.sim(a[[i]], b[[j]])
  r
}

pdf("git/sort_paper/FACS/spencerComparison/geneEnrichComparison.pdf",
  width=9, height=3.8)
par(mfrow=c(1,3))

a = overlap.matrix(bool.matrix.to.list(m.unnormalized > 0),
  bool.matrix.to.list(spencer.m))[,13:1]
plot.heatmap(a, main="Overlap of cells",
  cluster=FALSE, z.lim=c(0, max(a)/2))

b = overlap.matrix(facs.1, spencer.cluster)[,13:1]
plot.heatmap(b, main="Overlap of enriched genes",
  cluster=FALSE, z.lim=c(0, max(b)/2))

plot(as.vector(a), as.vector(b),
  xlab="Overlap of cells (Jaccard index)",
  ylab="Overlap of enriched genes (Jaccard index)", 
  pch=20, col="#00000080")

dev.off()

