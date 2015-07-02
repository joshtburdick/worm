/* A cluster browser. */

var clusterBrowser = {};

/* How each row is labelled. */
clusterBrowser.rowLabel = [
  {name:"gene", width:40},
  {name:"description", width:100},
  {name:"cluster", width:20}
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
    clusterBrowser.rowLabelSVG[i] = new SvgLabel(100, 20, 200);
    var g = clusterBrowser.rowLabelSVG[i].g;
    g.setAttributeNS(null, "transform", "translate(" + x + ",0)");
    h.appendChild(g);
    x = x + clusterBrowser.rowLabel[i].width;
  }

  // then, add column labels across the top
  columnLabels = new SvgLabel(a.arrayName.length, cellSize, 200);
  var g = columnLabels.g;
  for(var i=0; i<a.arrayName.length; i++)
    columnLabels.setText(i, a.arrayName[i]);
  g.setAttributeNS(null, "transform",
    "rotate(270," + x + ",0) translate(" + (x-200) + ",0)");
  document.getElementById("columnLabels").appendChild(g);

  // lastly, add the heatmap
  clusterBrowser.heatmapSVG = new SvgHeatmap(100, a.arrayName.length, cellSize);
  var g = clusterBrowser.heatmapSVG.g;
  g.setAttributeNS(null, "transform", "translate(" + x + ",0)");
  h.appendChild(g);
}

/* Updates the heatmap based on the page hash location. */
function updateHeatmap() {


}

/* If update is clicked, try to parse the gene in the gene field,
  and recenter on it. */
function updateClicked() {



}





