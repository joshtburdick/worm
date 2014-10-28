# Writes a large TSV file, listing clusters and TF-related enrichments.

library("hwriter")

source("git/utils.r")

clustering.dir = "git/cluster/hierarchical/"
output.dir = "git/sort_paper/plot/web/clusters"

# which clustering to write out
clustering.name = "hier.300.clusters"

# get all enriched motif and ChIP signals (not just significant cases)
motif.enriched.all = read.tsv(gzfile(paste0("git/sort_paper/tf/summary/motif/",
  clustering.name, "_merged.tsv.gz")))

chip.enriched.all = read.tsv(gzfile(paste0("git/sort_paper/tf/summary/chip/",
  clustering.name, "_unfiltered.tsv.gz")))
colnames(chip.enriched.all)[1] = "experiment"
chip.enriched.all$factor = sub("_.*$", "", chip.enriched.all$experiment)

# make ChIP factor names lower case (to help match them with
# motif names), and otherwise tweak them
i = tolower(chip.enriched.all$factor) %in% rownames(rr)
chip.enriched.all$factor[i] = tolower( chip.enriched.all$factor[i] )
chip.enriched.all$factor = sub("LIN-15B", "lin-15", chip.enriched.all$factor)

# Summary of per-cluster information.
cluster.tf.table = function() {

  # motif.enriched$group %in% c(1,2,3,52,286) 
  r = motif.enriched.all[ , c("group", "motif", "enrich", "p.corr") ]
# r = motif.enriched.all[ motif.enriched$group %in% as.character(c(1,2,3,52,286))  , c("group", "motif", "enrich", "p.corr") ]
  colnames(r) = c("group", "motif.id", "motif.enrich", "motif.p")

  m1 = unique(motif.ortholog[ , c("gene", "motif.id", "motif.name") ])
  r = merge(r, m1)

  # merge in cases in which a ChIP signal was enriched
  ch1 = chip.enriched.all[ , c("group","experiment","enrich","p.corr","factor") ]
  colnames(ch1) = c("group", "chip.experiment", "chip.enrich", "chip.p", "gene")
  r = merge(r, ch1, all.x=TRUE, all.y=TRUE)

  # add in correlation (when known)
  r$tf.corr = NA
  i = r$gene %in% colnames(tf.cluster.cor)
  # XXX might be cleaner to just convert tf.cluster.cor to a data.frame
  r[ i, "tf.corr" ] = tf.cluster.cor[
    cbind(match(r[i,"group"], rownames(tf.cluster.cor)),
      match(r[i,"gene"], colnames(tf.cluster.cor))) ]

  # variance in each cluster's enrichments
  cluster.enrich.var = apply(as.matrix(cluster.means), 1, var)
  r$cluster.enrich.var = cluster.enrich.var[ r$group ]

  r = r[ , c(1,2,11,10,6,4,5,7,8,9) ]
  colnames(r) = c("Cluster", "TF", "Cluster enrich var.",
    "TF-cluster corr.", "Motif", "Motif enrich", "Motif p",
    "ChIP", "ChIP enrich", "ChIP p")
  r$order = 5 * rank(-r$"Cluster enrich var.") +
    rank( - abs(r$"TF-cluster corr.") ) +
    rank(r$"Motif p") +
    rank(r$"ChIP p")

  r = r[ order(r$order) , ]

  r$"Cluster enrich var." = signif(r$"Cluster enrich var.", 3)
  r$"TF-cluster corr." = round(r$"TF-cluster corr.", 2)
  r$"Motif enrich" = round(r$"Motif enrich", 2)
  r$"Motif p" = signif(r$"Motif p", 2)
  r$"ChIP enrich" = round(r$"ChIP enrich", 2)
  r$"ChIP p" = signif(r$"ChIP p", 2)
  rownames(r) = NULL

  r
}

# Writes an index page.
# (FIXME For now, this will be a .tsv file.)
write.index = function() {

  b = ""

  s = paste0("<!DOCTYPE html>",
    "<html><head><title>Clustered FACS-sorted cells</title>",
    "<style type=\"text/css\">", css, "</style>", "</head>",
    "<body>", n, "</body></html>")
  

  cat(s, file=paste0(output.dir, "/", clustering.name, "/index.html"))
}

r = cluster.tf.table()
write.tsv(r, "git/sort_paper/network/clusterTF.tsv")

