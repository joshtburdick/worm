/* A cluster browser. */

var clusterBrowser = {};

clusterBrowser.numRows = 100;

// indices of gnees to show
clusterBrowser.p = [];

// for testing: initialize this arbitrarily
for(var i=0; i<100; i++)
  clusterBrowser.p[i] = i + 1234;

/* How each row is labelled. */
clusterBrowser.rowLabel = [
  {name:"gene", width:100},
  {name:"description", width:300}
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
console.log("added " + i + ", x = " + x);
    clusterBrowser.rowLabelSVG[i] = new SvgLabel(100, cellSize, 400);
    var g = clusterBrowser.rowLabelSVG[i].g;
    g.setAttributeNS(null, "transform", "translate(" + x + ",0)");
    h.appendChild(g);
    x = x + clusterBrowser.rowLabel[i].width;
  }

  // then, add column labels across the top
  columnLabelSVG = new SvgLabel(a.arrayName.length, cellSize, 200);
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

  // draw the initial set of genes
  redraw();

  // standardize (mean-center) the data
  xStandardized = a.x.map( vectorMath.standardize );
}

/* Updates the heatmap based on the page hash location. */
function redraw() {

  // update the heatmap
  h = clusterBrowser.heatmapSVG;
  p = clusterBrowser.p;
  for(var i=0; i<clusterBrowser.numRows; i++) {
    h.setRow(i, a.x[ p[i] ]);
  }

  // update row labels
  // XXX this shouldn't be hard-coded
  var nameLabel = clusterBrowser.rowLabelSVG[0];
  var descrLabel = clusterBrowser.rowLabelSVG[1];

  for(var i=0; i<clusterBrowser.numRows; i++) {
    nameLabel.setText(i, a.geneName[p[i]]);
    descrLabel.setText(i, a.descr[p[i]]);
  }


}

// Centers the clustering on a given gene (or genes; for now, just one gene)
// Args:
//   genes - a gene name, or a string of gene names separated by spaces
function setClustering(genes) {



  geneField.value = genes;

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
  goog.array.sort(p, function(i,j) {
    return r[j] - r[i];
  } );
console.log("after sort: p[0] = " + p[0]);

  // actually draw stuff, in that order
  drawClusters(graphics, a, p);

  // assuming all of the above worked, update the hash,
  // so that this is bookmarkable
//  location.hash = geneName1.toString();

//  graphics.render(canvas);
}

/* If update is clicked, try to parse the gene in the gene field,
  and recenter on it. */
function updateClicked() {



}





