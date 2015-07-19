# Writes an index page of all the clusters.

library("hwriter")

source("git/utils.r")
source("git/data/wormbase/wb.cluster.name.r")

clustering.dir = "git/cluster/hierarchical/"
output.dir = "git/sort_paper/plot/web/clusters"

tag = function(tag.name, x)
  paste0("<", tag.name, ">", x, "</", tag.name, ">")

# which clustering to write an index for
clustering.name = "hier.300.clusters"

# Utility to get the most significant enriched thing by cluster.
# Args:
#   e - a table of enrichments (with columns "group", "p.corr",
#     and column.name -- see below)
#   cl - the names of the clusters
#   column.name - the name of the column to return
#   pretty.print - function to apply to the most significantly
#     enriched thing
# Returns: vector of the same length as cl, containing how
#   many of the given things were enriched, and a description.
most.enriched = function(e, cl, column.name="group.name",
    pretty.print = function(a) a) {

  enriched.count = tapply(e[,column.name], e[,"cluster"], length)
  e = e[ order(e$p.corr) , ]
  e = e[ !duplicated(e$cluster) , ]
  e = tapply(e[,column.name], e[,"cluster"], function(x) x[1])

  r = rep(NA, length(cl))
  names(r) = cl
  r[ names(e) ] = paste0(enriched.count[names(e)], ": ", 
    sapply(e, pretty.print), "...")

  r[ is.na(r) ] = " "
  r
}

# Formats a cluster name somewhat.
cluster.name.format.old = function(a) {
  if (grepl("WBPaper00025032|WBPaper00036286", a))
    return(a)
  a = sub("^WBPaper[0-9]+:", "", a)
  a
}

# Formats a motif name.
motif.name.format = function(m) {
  # FIXME: these motif names are mostly lame
  paste0(motif.info[m,"motif.name"], " <img src=\"motifSvg/", m, ".svg\" ",
    "title=\"Logo for motif ", m, "\" width=\"107px\" height=\"24px\">")
}

# Creates an index table.
index.table = function() {

  a = cluster.means.1

  gene.count = as.vector(table(x$Cluster)[ rownames(a) ])

  # "colorize" the actual data
  for(j in 1:38) {
    a[,j] = sapply(a[,j], td.color)
  }

  # omit "Tissue"
  ao.enriched = ao.enriched[ ao.enriched$group.name != "Tissue" , ]
  colnames(hyperg.motif)[2] = "cluster"  # need to rename this

  motif.html = most.enriched(hyperg.motif, rownames(a),
    column.name = "motif", pretty.print = motif.name.format)

#  motif.html = sapply(motifs,
#    function(m) paste0("<img src=\"motifSvg/", m, ".svg\" ",
#      "title=\"Logo for motif ", m, "\" width=\"107px\" height=\"24px\">"))

  colnames(hyperg.chip)[2] = "cluster"   # again, renaming this
  chip = most.enriched(hyperg.chip, rownames(a), column.name="experiment")
  chip = sub("-post-early-Embryos", "", chip)   # XXX

  # tack on annotation
  cluster.links = html.link(rownames(a),
    sapply(rownames(a),
      function(cl) paste0(clustering.name, "/", cl, ".html"))) 
  a = cbind(a, data.frame(
    "Cluster"=tag("td", cluster.links),
    "Number of genes"=tag("td", gene.count),
    "Mean expr. (RPM)"=tag("td", round(mean.max.rpm[ rownames(a) ], 1)),
    "Anatomy enriched"=tag("td", most.enriched(ao.enriched, rownames(a))),
#    "Phenotype enriched"=tag("td", most.enriched(phenotype.enriched, rownames(x))),
    "Other expression cluster"=tag("td", most.enriched(cluster.enriched, rownames(a), column.name="group", pretty.print=cluster.name.format)),
    "Motif"=tag("td", motif.html),
    "ChIP"=tag("td", chip),
  check.names=FALSE))

  # column names (some rotated using CSS hacking)
  col.names = colnames(a)
  col.names[1:41] = sapply(col.names[1:41], rotate.text)

  table.header = paste(c("<thead><tr>",
    sapply(col.names, function(a) paste0("<th width=\"20px\" height=\"200px\">", a, "</th>")),
    "</thead></tr>\n"), collapse="")

  # we assume contents of x already have been formatted
  table.rows.1 = apply(a, 1, function(r)
    paste0("<tr height=12px>", paste(r, collapse=""), "</tr>\n", collapse=""))
  table.rows = paste0("<tbody class=expr-tbody>", paste(table.rows.1, collapse=""), "</tbody>")

  paste0("<table class=expr-table>", table.header, table.rows, "</table>")
}

# Writes out an index page.
write.index.page = function() {
  # hack to read in a text file as a string
  css = paste0(read.table("git/sort_paper/plot/web/clusterCss.txt", sep="~")[,1], collapse="")

  h = index.table()

  s = paste0("<!DOCTYPE html>",
    "<html><head><title>", "Cluster index</title>",
    "<style type=\"text/css\">", css, "</style>",
    "<script>", "</script>", "</head>",
    "<body>",
    "<h1>Cluster index</h1>",
    h, "</body></html>")
  cat(s, file=paste0(output.dir, "/", 
    "index.html"))
}

write.index.page()

