/* A cluster browser.
TODO:
20150705
- After going to Wormbase, coming back requires pressing
  "back" twice. Presumably linking by Wormbase ID would fix this.
- Presumably, the heatmap should be scrolled to the top after
  jumping to a gene (otherwise, it's confusing, because the gene
  you jumped to isn't visible.)
- should allow starting on an arbitrary gene (e.g. pha-4).
*/

var clusterBrowser = {};

/* how many genes to show */
clusterBrowser.numRows = 100;

/* indices of genes to show */
clusterBrowser.p = [];

/* initialize this to "unsorted" */
for(var i=0; i<a.x.length; i++)
  clusterBrowser.p[i] = i;

/* How each row is labelled. */
clusterBrowser.rowLabel = [
  {name:"gene", width:100, colors:["#fff", "#eee"]},
  {name:"description", width:300, colors:["#eee","#ddd"]},
  {name:"cluster", width:40, colors:["#fff","#eee"]}
];


/* Initializes the page: adds the relevant CSS. */
function initPage() {
  var cellSize = 20;

  var h = document.getElementById("heatmap");

  // first, add row labels (from left to right, so that
  // they'll overlap, since I don't know how to clip these)
  clusterBrowser.rowLabelSVG = [];
  var x = 0;
  for(var i=0; i<clusterBrowser.rowLabel.length; i++) {
// console.log("added " + i + ", x = " + x);
    clusterBrowser.rowLabelSVG[i] =
      new SvgLabel(100, cellSize, 400, clusterBrowser.rowLabel[i].colors);
    var g = clusterBrowser.rowLabelSVG[i].g;
    g.setAttributeNS(null, "transform", "translate(" + x + ",0)");
    h.appendChild(g);
    x = x + clusterBrowser.rowLabel[i].width;
  }

  // then, add column labels across the top
  columnLabelSVG = new SvgLabel(a.arrayName.length, cellSize, 200,
    ["#fff","#eee"]);
  var g = columnLabelSVG.g;
  for(var i=0; i<a.arrayName.length; i++)
    columnLabelSVG.setText(i, a.arrayName[i]);
  g.setAttributeNS(null, "transform",
    "rotate(270," + x + ",0) translate(" + (x-200) + ",0)");
  document.getElementById("columnLabels").appendChild(g);

  // lastly, add the heatmap
  clusterBrowser.heatmapSVG = new SvgHeatmap(100, a.arrayName.length, cellSize);
  var g = clusterBrowser.heatmapSVG.g;
  g.setAttributeNS(null, "transform", "translate(" + x + ",0)");
  h.appendChild(g);

  // standardize (mean-center) the data
  xStandardized = a.x.map( vectorMath.standardize );

  // if there isn't a gene defined, go to pha-4
  if (location.hash.length <= 1) {
    location.hash = "pha-4";
  }

  // draw the initial set of genes
  redraw();

  window.onhashchange = setClustering;
}

/* Updates the heatmap based on the page hash location. */
function redraw() {

  // update the heatmap
  h = clusterBrowser.heatmapSVG;
  p = clusterBrowser.p;
  for(var i=0; i<clusterBrowser.numRows; i++) {
    h.setRow(i, a.x[ p[i] ]);
  }

  // update row labels (and links)
  // XXX this shouldn't be hard-coded
  var nameLabel = clusterBrowser.rowLabelSVG[0];
  var descrLabel = clusterBrowser.rowLabelSVG[1];
  var clusterLabel = clusterBrowser.rowLabelSVG[2];

  for(var i=0; i<clusterBrowser.numRows; i++) {
    nameLabel.setText(i, a.geneName[p[i]]);
    nameLabel.setLink(i,
      "#" + a.geneName[p[i]],
      "recenter on " + a.geneName[p[i]]);

    descrLabel.setText(i, a.descr[p[i]]);
    descrLabel.setLink(i,
      "http://www.wormbase.org/db/get?name=" + a.geneName[p[i]] + ";class=Gene",
      "Wormbase on " + a.geneName[p[i]]);

    clusterLabel.setText(i, a.cluster[p[i]]);
    clusterLabel.setLink(i,
      "../clusters/hier.300.clusters/" + a.cluster[i] + ".html",
      "go to cluster " + a.cluster[i]);
  }
}

// Centers the clustering on a given gene (or genes; for now, just one gene)
function setClustering() {

  // remove the initial #
  genes = location.hash.substring(1);

  // look up gene names
  geneIndex = new Array();
  geneName1 = new Array();
  s = genes.split(" ");
  for(i=0; i<s.length; i++) {
    if (a.geneNameToIndex[ s[i] ] != undefined) {
      geneIndex.push( a.geneNameToIndex[s [i] ]);
      geneName1.push( s[i] );
    }
  }

  // if we didn't find a gene, return
  if (geneIndex.length == 0) {
    return;
  }

  // update the field
  document.getElementById("genes").value = genes;

  // FIXME: show an error if no gene names are legit
  // compute center of all selected genes
  // XXX for now, just using the first gene

  xCenter = xStandardized[ geneIndex[0] ];

  console.log("geneIndex = " + geneIndex);

  // compute correlations with this
  r = xStandardized.map( function(x) { return vectorMath.dot(x, xCenter); } );
console.log("updated r");

console.log("before sort: p[0] = " + p[0]);

  // update the ordering, p
  // FIXME force the selected genes to be first (although presumably
  // they'll already be near the start)
  p = p.sort(function(i,j) {
    return r[j] - r[i];
  } );
console.log("after sort: p[0] = " + p[0]);

  // redraw the genes, in that order
  redraw();
}

/* If update is clicked, try to parse the gene in the gene field,
  and recenter on it. */
function updateClicked() {
  console.log("update clicked: field is " +
    document.getElementById("genes").value);
  location.hash = document.getElementById("genes").value;
}





