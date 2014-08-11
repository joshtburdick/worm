# Writes out clusters of data.

library("hwriter")

source("git/utils.r")

clustering.dir = "git/cluster/hierarchical/"
output.dir = "git/sort_paper/plot/web/clusters"

# which clustering to write out
clustering.name = "hier.300.clusters"

if (FALSE) {
go.enriched = read.table(
  "git/sort_paper/enrichment/geneOntology/hier.300.clusters.tsv",
  sep="\t", quote="", header=TRUE, as.is=TRUE)
go.enriched = go.enriched[ go.enriched$Count >= 2 , ]
}
system(paste0("mkdir -p ", output.dir, "/", clustering.name))

# read in data to cluster
x = read.table(paste0(clustering.dir, "/", clustering.name, "/cluster.cdt"),
  sep="\t", quote="", fill=TRUE, header=TRUE, check.names=FALSE, as.is=TRUE)
x = x[ c(-1,-2) , ]
for(j in c(12:51)) {
  x[ , j ] = as.numeric( x[ , j ])
}

# enriched anatomy terms and other expression clusters
ao.enriched = read.tsv(paste0(
  "git/sort_paper/enrichment/summary/anatomyEnrichment/",
  clustering.name, ".tsv"))
cluster.enriched = read.tsv(paste0(
  "git/sort_paper/enrichment/summary/wormbaseCluster/",
  clustering.name, ".tsv"))

# enriched motif and ChIP signals
motif.enriched = read.tsv(paste0("git/sort_paper/tf/summary/motif/",
  clustering.name, ".tsv"))
chip.enriched = read.tsv(paste0("git/sort_paper/tf/summary/chip/",
  clustering.name, ".tsv"))



# Utility to convert numbers to colors.
# Args:
#   neg.hue, pos.hue - the negative and positive hues
#   max.val - the value corresponding to "fully saturated"
#   x - the actual number
number.to.color.hsv = function(neg.hue, pos.hue, max.val, x) {
  x = x / max.val
  x[ x > 1 ] = 1
  x[ x < -1 ] = -1

  hsv(ifelse(x < 0, neg.hue, pos.hue), 1, abs(x))
}

# HTML for an empty table cell, colored some color based on
# a number (which shows up on mouseover).
td.color = function(x) {
  # FIXME use HTML5
  paste0("<td title=\"", signif(x,3), "\" bgcolor=\"",
    number.to.color.hsv(1/6, 4/6, 2, x),
    "\" > </td>\n")
}

# HTML5 to rotate text 90 degrees CCW (assuming the above
# stylesheet has been loaded.)
rotate.text = function(s) {
  paste0("<div class=\"vertical-text\">", s, "</div>", collapse="")


}

# HTML for a table of enriched anatomy terms / clusters.
enriched.table = function(name, a) {
  if (nrow(a) == 0)
    return("")

  a = a[ , c(2,3,4,6,7,8) ]
  a = a[ order(a$p.corr) , ]
  a[ , c(5,6) ] = signif(a[ , c(5,6) ], 2)

  # XXX ideally hwrite should write the first line in <th>
  colnames(a) = sapply(c("Group", "Group name", "Group size",
    "Number in cluster", "p", "FDR corrected p"),
    function(a) paste0("<b>", a, "</b>"))

  paste0("<h3>", name, "</h3>", hwrite(a, row.names=FALSE))
}

# HTML for a table of enriched GO terms.
go.enriched.table = function(a) {
  a = a[ order(a$p.corr) , ]
  results = data.frame("GO term" = a$Term, "GO ID" = a$GOID,
    "Number of genes in cluster" = paste0(a$Count, "/", a$Size),
    "FDR-corrected p-value" = signif(a$p.corr, 2), check.names=FALSE)
  hwrite(results)
}

# HTML for a table of enriched motifs / ChIP terms.
motif.chip.table = function(a, motif.or.chip) {
  a = a[ order(a$p.corr) , ]

  # XXX ideally hwrite should write the first line in <th>
  colnames(a) = sapply(
    c(ifelse(motif.or.chip=="motif", "Motif", "ChIP experiment"),
    "Group", "Upstream distance (kb)", "Conservation",
    "Motif score", "Motifs in cluster", "Motifs in background",
    "Sequence in cluster (bp)",
    "Sequence in background (bp)", "Enrichment", "Chi-squared value",
    "p", "FDR corrected p"),
    function(a) paste0("<b>", a, "</b>"))

  a[ , c(10:11) ] = round(a[,c(10:11)], 1)
  a[ , c(12:13) ] = signif(a[,c(12:13)], 2)

  if (motif.or.chip == "chip") {
    a = a[ , -5 ]
  }
  a = a[ , -2 ]

  hwrite(a, row.names=FALSE)
}

# HTML for a formatted table of expression numbers.
expr.table = function(x, cl) {

  # FIXME this should include links
  for(j in 1:11) {
    x[,j] = sapply(x[,j], function(a) paste0("<td>", a, "</td>"))
  }

  # "colorize" the actual data
  for(j in 12:51) {
    x[,j] = sapply(x[,j], td.color)
  }
  x = x[ , c(12:51, 3, 4) ]

  col.names = colnames(x)
  col.names[1:40] = sapply(col.names[1:40], rotate.text)

  table.header = paste(c("<thead><tr>",
    sapply(col.names, function(a) paste0("<th width=\"20px\" height=\"200px\">", a, "</th>")),
    "</thead></tr>\n"), collapse="")

  # we assume contents of x already have been formatted
  table.rows.1 = apply(x, 1, function(r)
    paste0("<tr height=12px>", paste(r, collapse=""), "</tr>\n", collapse=""))
  table.rows = paste0("<tbody class=expr-tbody>", paste(table.rows.1, collapse=""), "</tbody>")

  paste0("<table class=expr-table>", table.header, table.rows, "</table>")
}

# Writes out one cluster.
write.cluster = function(x, cl) {
  write.status(cl)

  x1 = x[ x$Cluster == cl , ]
#  x1 = x[1:10,]

  # hack to read in a text file as a string
  css = paste0(read.table("git/sort_paper/plot/web/clusterCss.txt", sep="~")[,1], collapse="")
#  css = ""  # XXX hack to disable CSS

  h = paste0(
    "<h2>Cluster ", cl, "</h2>",
    "<h3>Anatomy terms</h3>",
      enriched.table("Anatomy terms enriched",
        ao.enriched[ ao.enriched$cluster==cl , ]),
    "<h3>GO terms enriched (FIXME)</h3>",
      enriched.table("Expression clusters enriched",
        cluster.enriched[ cluster.enriched$cluster==cl , ]),

    "<h3>Motifs enriched</h3>",
    motif.chip.table(motif.enriched[motif.enriched$group == cl,], "motif"),
    "<h3>ChIP peaks enriched</h3>",
    motif.chip.table(chip.enriched[chip.enriched$group == cl,], "chip"),
    "<h3>Expression clusters</h3>", expr.table(x1)
  )

  s = paste0("<!DOCTYPE html>",
    "<html><head><title>", "Cluster ", cl, "</title>",
    "<style type=\"text/css\">", css, "</style>", "</head>",
    "<body>", h, "</body></html>")
  cat(s, file=paste0(output.dir, "/", clustering.name, "/",
    cl, ".html"))
}

for(cl in sort(unique(x$Cluster))) {
  write.cluster(x, cl)
}
