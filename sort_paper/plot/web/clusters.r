# Writes out clusters of data.

library("hwriter")

source("git/utils.r")
source("git/tf/motif/motifName.r")
source("git/plot/web/html.r")

clustering.dir = "git/cluster/hierarchical/"
output.dir = "git/sort_paper/plot/web/clusters"

# which clustering to write out
clustering.name = "hier.300.clusters"

# WormBase gene IDs (for linking to there)
wb.gene.id = read.csv(gzfile("data/wormbase/geneIDs.WS224.csv.gz"),
  header=FALSE, as.is=TRUE)
colnames(wb.gene.id) = c("wb.gene", "name", "clone.id")

go.enriched = read.table(
  "git/sort_paper/enrichment/geneOntology/hier.300.clusters.tsv",
  sep="\t", quote="", header=TRUE, as.is=TRUE)
go.enriched = go.enriched[ go.enriched$Count >= 2 , ]

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
phenotype.enriched = read.tsv(paste0(
  "git/sort_paper/enrichment/summary/phenotype/",
  clustering.name, ".tsv"))

# enriched motif and ChIP signals
motif.enriched = read.tsv(paste0("git/sort_paper/tf/summary/motif/",
  clustering.name, ".tsv"))
chip.enriched = read.tsv(paste0("git/sort_paper/tf/summary/chip/",
  clustering.name, ".tsv"))
colnames(chip.enriched)[1] = "experiment"
chip.enriched$factor = sub("_.*$", "", chip.enriched$experiment)

# information about orthologs
motif.ortholog = read.table("git/tf/motif.ortholog.2.tsv",
  as.is=TRUE, header=TRUE)
motif.filter = read.tsv("git/tf/motif/motifFilter.tsv")
# add representative motif name, from the motif clustering
motif.ortholog$canonical.motif =
  motif.filter[ motif.ortholog$motif, "canonical.name" ]
motif.ortholog = motif.ortholog[ !is.na(motif.ortholog$canonical.motif) , ]
# for each motif, list of potential orthologs
orthologs.by.motif = by(motif.ortholog$gene, motif.ortholog$canonical.motif,
  function(x) {
    x = as.character(x)
    nhrs = grep("nhr", x)
    if (length(x) >= 5 && length(nhrs) >= 5) {
      return(unique(c(grep("nhr", x, invert=TRUE, value=TRUE), paste(length(nhrs), "NHRs"))))
    }

    return(unique(as.character(x))) 
  }
)

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
    number.to.color.hsv(194/360, 60/360, 2, x),
    "\" > </td>\n")
}

# HTML5 to rotate text 90 degrees CCW (assuming the above
# stylesheet has been loaded.)
rotate.text = function(s) {
  paste0("<div class=\"vertical-text\">", s, "</div>", collapse="")

}

# HTML for a table of enriched anatomy terms / clusters.
enriched.table = function(a, link.f=NULL) {
  if (nrow(a) == 0)
    return("none found")

#  a = a[ , c(2,3,4,6,7,8) ]

  if (!is.null(link.f)) {
    a[,"group.name"] =
      html.link(a[,"group.name"], link.f(a[,"group"]))
  }
  a = a[ , c("group.name", "num.intersect", "enrichment", "p.corr") ]
  a = a[ order(a$p.corr) , ]

#  a[ , c(5,6) ] = signif(a[ , c(5,6) ], 2)

  # XXX ideally hwrite should write the first line in <th>
  colnames(a) = sapply(c("Group name",
    "Number in cluster", "Enrichment", "FDR corrected p"),
    function(a) paste0("<b>", a, "</b>"))

  hwrite(a, row.names=FALSE)
}

# HTML for a table of enriched GO terms.
go.enriched.table = function(a) {
  if (nrow(a) == 0)
    return("none found")

  a = a[ order(a$p.corr) , ]
  results = data.frame(   #"GO term" = a$Term, "GO ID" = a$GOID,
    "GO term" = html.link(a$Term,
      paste0("http://www.wormbase.org/db/get?name=",
        a$GOID, ";class=go_term")),
    "Number of genes" = a$Count,
#    "Number of genes in cluster" = paste0(a$Count, "/", a$Size),
    "FDR-corrected p-value" = signif(a$p.corr, 2), check.names=FALSE)
  colnames(results) = sapply(colnames(results),
    function(a) paste0("<b>", a, "</b>"))
  hwrite(results)
}

# HTML for a table of enriched motifs.
motif.table = function(a) {
  if (nrow(a) == 0)
    return("none found")

  a = a[ order(a$p.corr) , ]

  a$logo = sapply(a$motif,
    function(m) paste0("<img src=\"../motifSvg/", m, ".svg\" ",
      "alt=\"Logo for motif ", m, "\" width=\"120px\" height=\"27px\">"))
  a$enrich = round(a$enrich, 2)
  a$p.corr = signif(a$p.corr, 2)
  a$ortholog = sapply(a$motif,
    function(m) {
      g = orthologs.by.motif[[ m ]]
      # if there are many TFs, show the entire list as alt text
      if (length(g) > 21) {
        all.g = paste(g, collapse=" ")
        g = c(g[1:20], paste("and", length(g)-20, "others"))

        paste(paste(g, collapse=" "),
          paste0("<br><span title=\"", all.g, "\">[full list]</span>"))
      }
      else {
        paste(g, collapse=" ")
      }
    })

  a = a[ , c("motif", "logo", "ortholog", "motifs.cluster", "enrich", "p.corr") ]

  # use human-readable motif names
  a$motif = motif.name[ a$motif ]

  # XXX ideally hwrite should write the first line in <th>
  colnames(a) = sapply(
    c("Motif", "Logo", "Possible orthologs",
      "Number of motifs in cluster", "Enrichment", "FDR corrected p"),
    function(a) paste0("<b>", a, "</b>"))

  hwrite(a, row.names=FALSE)
}


# HTML for a table of ChIP terms.
chip.table = function(a) {
  if (nrow(a) == 0)
    return("none found")

  a = a[ order(a$p.corr) , ]

  a = a[ , c("factor", "experiment", "motifs.cluster", "enrich", "p.corr") ]

  a$enrich = round(a$enrich, 2)
  a$p.corr = signif(a$p.corr, 2)
  colnames(a) = sapply(
    c("Gene", "Experiment", "Number of upstream peaks", "Enrichment",
    "FDR corrected p"),
    function(a) paste0("<b>", a, "</b>"))

  hwrite(a, row.names=FALSE)
}

# HTML for a formatted table of expression numbers.
expr.table = function(x, cl) {

# add link out to WormBase
#  x[,3] = sapply(x[,3],
#    function(g) html.link(g,
#      paste0("http://www.wormbase.org/db/get?name=", g, ";class=Gene")))
  wb.ids = x[,3]
  i = x[,3] %in% wb.gene.id[,3]
  wb.ids[i] = wb.gene.id[match(wb.ids[i], wb.gene.id[,3]), 1]  
  i = x[,3] %in% wb.gene.id[,2]
  wb.ids[i] = wb.gene.id[match(wb.ids[i], wb.gene.id[,2]), 1]  

  x[,3] = html.link(x[,3],
    paste0("http://www.wormbase.org/db/get?name=", wb.ids, ";class=Gene"))

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

# JavaScript for "go to next page."
next.page.js = function(cl) {
  cl.next = cl + 1
  if (cl.next == 301) {
    cl.next = 1
  }

  paste0(
  "function keyHandler(event) {",
  " if (event.keyCode == 32) {",
  "   window.location = \"", cl.next, ".html\"; ",
  " }",
  "}")
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
    "<h3>Expression</h3>", expr.table(x1),
    "<h3>Anatomy terms enriched</h3>",
    enriched.table(ao.enriched[ ao.enriched$cluster==cl , ],
      function(ao)
        paste0("http://www.wormbase.org/db/get?name=",
        ao, ";class=Anatomy_term")),
    "<h3>Phenotypes enriched</h3>",
    enriched.table(phenotype.enriched[ phenotype.enriched$cluster==cl , ],
      function(a)
        paste0("http://www.wormbase.org/db/get?name=",
        a, ";class=Phenotype")),
    "<h3>GO terms enriched</h3>",
    go.enriched.table(go.enriched[ go.enriched$cluster==cl , ]),
    "<h3>Expression clusters enriched</h3>",
    enriched.table(cluster.enriched[ cluster.enriched$cluster==cl , ],
      function(ao)
        paste0("http://www.wormbase.org/db/get?name=",
        ao, ";class=Expression_cluster")),
    "<h3>Motifs enriched</h3>",
    motif.table(motif.enriched[motif.enriched$group == cl,]),
    "<h3>ChIP peaks enriched</h3>",
    chip.table(chip.enriched[chip.enriched$group == cl,])
  )

  s = paste0("<!DOCTYPE html>",
    "<html><head><title>", "Cluster ", cl, "</title>",
    "<style type=\"text/css\">", css, "</style>",
    "<script>", next.page.js(cl), "</script>", "</head>",
    "<body onkeyup=\"keyHandler(event)\">",
    h, "</body></html>")
  cat(s, file=paste0(output.dir, "/", clustering.name, "/",
    cl, ".html"))
}

# for(cl in c(1,2,3,52,286)) {
for(cl in sort(unique(x$Cluster))) {
  write.cluster(x, cl)
}

