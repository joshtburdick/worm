
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

// Utility to color rows alternating colors.
function alternatingFill(i) {
  if (i % 2 == 0)
    return new goog.graphics.SolidFill('#f8f8f8');
  else
    return new goog.graphics.SolidFill('#e8e8e8');
}

// Draws a column of labels for the heatmap.
function drawLabels(g, bounds, s) {
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
function drawHeatmap(g, bounds, a) {
  colorMapping = exprToColor(exprRange, [0,0,255], [0,0,0], [255,255,0]);

  var x0 = bounds.left;
  var y0 = bounds.top;
  for(i=0; i<numRows; i++) {
    for(j=0; j<a.arrayName.length; j++) {
      var f = colorMapping(a.x.getValueAt(p[i],j));
      g.drawRect(x0 + cellSize * j, y0 + cellSize * i,
        cellSize + 1, cellSize + 1,
        null, new goog.graphics.SolidFill(f));
    }
  }
}

// Draws the clustering.
function drawClusters(g, a) {

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
  drawLabels(g, labelBounds, a.geneName);

  // draw the heatmap
//  heatmapBounds = new goog.math.Rect(
//    bounds.left + 250, bounds.top + 100, bounds.width, bounds.height);
  drawHeatmap(g, heatmapBounds, a);

}

// event handler for when a gene is clicked
function imageMouseDown(e) {
  loc = new goog.math.Coordinate(e.offsetX, e.offsetY);

  if (!(labelBounds.contains(loc)))
    return;

  i = Math.trunc( (e.offsetY - heatmapBounds.top) / cellSize );

//  alert("mouse clicked on " + a.geneName[ p[i] ]);

  var r = calcCorrelations(p[i]);
  p = r.p;
  drawClusters(graphics, a);

}

// This is initially called when the page loads.
function showClusters() {
  a = readCDT("cluster.cdt");

  // ordering of the genes
  p = new Array(numRows);
  for(i=0; i<numRows; i++)
    p[i] = i;

  xStandardized = standardizeRows(a.x);

  graphics = goog.graphics.createGraphics('1200px', '2000px');
  canvas = goog.dom.$('heatmap');
  graphics.render(canvas);

  drawClusters(graphics, a);

  goog.events.listen(canvas, goog.events.EventType.MOUSEDOWN,
    imageMouseDown);

}



