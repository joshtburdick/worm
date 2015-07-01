/* A cluster browser. */

var clusterBrowser = {};



/* Initializes the page: adds the relevant CSS. */
function initPage() {
  var cellSize = 20;
  var svg = document.getElementById("columnLabels");

  var g = document.createElementNS(svgns, "g");
  g.setAttributeNS(null, "width", "1500");
  g.setAttributeNS(null, "height", "1500");
  g.setAttributeNS(null, "transform", "rotate(270,150,0) translate(0,0)")

  columnLabels = new SvgLabel(g, a.arrayName.length, cellSize);
  for(i=0; i<a.arrayName.length; i++)
    columnLabels.setText(i, a.arrayName[i]);
  svg.appendChild(g);
  console.log("added labels");

  var h = document.getElementById("heatmap");
  clusterBrowser.heatmap = new SvgHeatmap(h,
    100, a.arrayName.length, cellSize);



}

/* Updates the heatmap based on the page hash location. */
function updateHeatmap() {


}

/* If update is clicked, try to parse the gene in the gene field,
  and recenter on it. */
function updateClicked() {



}





