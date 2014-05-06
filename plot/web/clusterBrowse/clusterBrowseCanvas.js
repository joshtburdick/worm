// Draws browsable clusters.

goog.require('goog.Uri');
goog.require('goog.color');
goog.require('goog.dom');
goog.require('goog.events');
goog.require('goog.graphics');
goog.require('goog.math');

// Size of a heatmap cell.
cellSize = 20;

// Number of rows to show.
numRows = 200;

// bound for the expression range coloring
exprRange = 2;

// index of gene to "center" the clustering on
clusterCenterIndex = 13;

// Reads in clustering data, in CDT format.
function readCDT(url) {

  var s = readTSV(url);

  // number of genes
  var m = s.length - 3;
  var arrayName = s[0].slice(4);
  var geneName = new Array(m);
  var geneNameToIndex = {};
  var descr = new Array(m);
  var x = new Array(m);

  // loop through the lines, filling in gene info and expression data
  for(i=0; i<m; i++) {
    var f = s[i+3];
    geneName[i] = f[1];
    geneNameToIndex[f[1]] = i;
    descr[i] = f[2];
    x[i] = f.slice(4).map(parseFloat);
  }

  return { arrayName:arrayName, geneName:geneName,
    geneNameToIndex:geneNameToIndex, descr:descr, x:x };
}

// Draws a column of labels for the heatmap.
function drawLabels(g, bounds, s, p) {
  var labelFont = new goog.graphics.Font(15, 'Times');
  var fill = new goog.graphics.SolidFill('black');

  for(i=0; i<numRows; i++) {
    g.drawRect(bounds.left, bounds.top + cellSize * i,
      bounds.width, cellSize, null, alternatingFill(i));
    g.drawText(s[p[i]], bounds.left, bounds.top + cellSize * i,
      bounds.width, cellSize, 'left', 'center', labelFont, null, fill);
  }
}

// Draws the array labels (along the top).
function drawArrayLabels(g, bounds, arrayName) {
  var labelFont = new goog.graphics.Font(15, 'Times');
  var fill = new goog.graphics.SolidFill('black');

  for(j=0; j<arrayName.length; j++) {
    g.drawRect(bounds.left + cellSize * j, bounds.top,
      cellSize, bounds.height, null, alternatingFill(j));
    var x1 = bounds.left + cellSize * j + cellSize / 2;
    g.drawTextOnLine(arrayName[j],
      x1, bounds.top + bounds.height - 3, x1, bounds.top,
      'left', labelFont, null, fill);
  }
}

// Draws a heatmap of the clusters.
function drawHeatmap(g, bounds, a, p) {

  // blue-black-yellow color gradient
  colorMapping = exprToColor(exprRange, [0,0,255], [0,0,0], [255,255,0]);

  var x0 = bounds.left;
  var y0 = bounds.top;
  for(i=0; i<numRows; i++) {
    for(j=0; j<a.arrayName.length; j++) {
      var f = colorMapping(a.x[ p[i] ][ j ]);
      g.drawRect(x0 + cellSize * j, y0 + cellSize * i,
        cellSize + 1, cellSize + 1,
        null, new goog.graphics.SolidFill(f));
    }
  }
}

// Draws the clustering.
// Args:
//   g - graphics context on which to draw
//   a - object containing the data to render
//   p - permutation of rows
function drawClusters(g, a, p) {

  console.log("in drawClusters: p[0] =" + p[0]);
  bounds = new goog.math.Rect(0, 0,
    1000, 2000);  // FIXME shouldn't be hardcoded

  // where to put the actual heatmap
  heatmapBounds = new goog.math.Rect(250, 150, 950, 1800);

  // draw various labels
  drawArrayLabels(g, new goog.math.Rect(
    heatmapBounds.left, bounds.top, bounds.width, heatmapBounds.top),
    a.arrayName);

  labelBounds = new goog.math.Rect(
    bounds.left, heatmapBounds.top, heatmapBounds.left, numRows * cellSize);
  drawLabels(g, labelBounds, a.geneName, p);

  // draw the heatmap
//  heatmapBounds = new goog.math.Rect(
//    bounds.left + 250, bounds.top + 100, bounds.width, bounds.height);
  drawHeatmap(g, heatmapBounds, a, p);
}


// Centers the clustering on a given gene (or genes; for now, just one gene)
// Args:
//   genes - a gene name, or a string of gene names separated by
//     commas or spaces
function setClustering(genes) {
  geneField.value = genes;

  console.log("genes = " + genes);

  // look up gene names
  geneIndex = new Array();
  geneName1 = new Array();
  s = genes.split(",");
  for(i=0; i<s.length; i++) {
    if (a.geneNameToIndex[ s[i] ] != undefined) {
      geneIndex.push( a.geneNameToIndex[s [i] ]);
      geneName1.push( s[i] );
    }
  }

  // FIXME: check for legit gene names


  // center of all selected genes
  // FIXME: for now, just uses the first one
  xCenter = xStandardized[ geneIndex[0] ];

console.log("geneIndex = " + geneIndex);

  // compute correlations with this
  r = xStandardized.map( function(x) { return vectorMath.dot(x, xCenter); } );
console.log("updated r");

console.log("before sort: p[0] = " + p[0]);
  // update the ordering, p
  goog.array.sort(p, function(i,j) { return r[j] - r[i]; } );
console.log("after sort: p[0] = " + p[0]);

  // actually draw stuff, in that order
  drawClusters(graphics, a, p);

//  graphics.render(canvas);

  // assuming all of the above worked, update the hash,
  // so that this is bookmarkable
//  location.hash = geneName1.toString();
}

// event handler for when a gene is clicked
function imageMouseUp(e) {
  // first, get mouse click location relative to this canvas
  // see, e.g.
  // http://stackoverflow.com/questions/55677/how-do-i-get-the-coordinates-of-a-mouse-click-on-a-canvas-element
  canvasLoc = document.getElementById("heatmap").getBoundingClientRect();
console.log("canvas offset = " + canvasLoc.left + " " + canvasLoc.top);
  loc = new goog.math.Coordinate(e.clientX - canvasLoc.left, e.clientY - canvasLoc.top);
console.log("loc = " + loc.toString());
  if (!(labelBounds.contains(loc)))
    return;

  i = Math.trunc( (loc.y - heatmapBounds.top) / cellSize );
  if (i<0) i=0;

  console.log("mouse clicked on " + a.geneName[ p[i] ]);

//  setClustering( a.geneName[ p[i] ] );


  uri = new goog.Uri(location.search);
  uri.setParameterValue("gene", a.geneName[ p[i] ] );
  location.search = uri;
}

// This is initially called when the page loads.
function clusterBrowseInit() {
  geneField = document.getElementById("genes");
  geneField.focus();

console.log("in clusterBrowseInit()");
  // set up for drawing
  graphics = goog.graphics.createGraphics('1200px', '2000px');
  canvas = goog.dom.$('heatmap');

  uri = new goog.Uri(location.search);
  genes = uri.getParameterValue("gene");

  // update the gene names in the text box
  geneField.value = genes;

  // read in the data, and standardize (mean-center) it
  a = readCDT("cluster.cdt");
  xStandardized = a.x.map( vectorMath.standardize );

console.log("standardized rows");

  minExpr = uri.getParameterValue("minExpr");
  if (minExpr == undefined) {
    minExpr = 1;
  }

/*
  // also get gene(s) to center on; I'm using the hash
  // portion to avoid reloading the page all the time.
// XXX deprecated
  genes =  location.hash;
  if (location.hash.length >= 2) {
    genes = location.hash.slice(1);
  }
  else {
    genes = "pha-4";
  }
*/

  // ordering of the genes
  p = new Array(a.geneName.length);
  for(i=0; i<a.geneName.length; i++)
    p[i] = i;

  setClustering(genes);
//  drawClusters(graphics, a);

  // request mouse events (now that the ordering p is defined)
  goog.events.listen(canvas, goog.events.EventType.MOUSEUP,
    imageMouseUp);

/*
  window.location.watch(
    'hash',
    function(id,oldVal,newVal){
      console.log("the window's hash value has changed from "+oldval+" to "+newVal);
//      setClustering(newVal.slice(1) );
    }
);
*/

  graphics.render(canvas);
}

// This should be called when the form's "update" button is clicked.
function updateClicked() {

  console.log("update clicked: genes = " + geneField.value);
//  setClustering(geneField.value);

  // update the gene names in the text box; this, in turn,
  // will trigger going to the appropriate page
  uri = new goog.Uri(location.search);
  uri.setParameterValue("gene", geneField.value);
  location.search = uri;
}

