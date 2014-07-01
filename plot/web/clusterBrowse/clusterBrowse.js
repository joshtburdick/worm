
goog.require('goog.color');
goog.require('goog.dom');
goog.require('goog.graphics');
goog.require('goog.math.Matrix');

// Renders some rows of the matrix as an HTML table.
// Args:
//   a - the data, as an object with fields
//     arrayName, geneName, descr, and x (as
//     returned by readCDT())
//   r - integer indices of which rows to show
// Returns: an HTML table which displays that.
function showClustersAsTable(a, r) {

  loColor = [0,0,255];
  zeroColor = [0,0,0];
  hiColor = [255,255,0];

  colorMapping = exprToColor(1, loColor, zeroColor, hiColor);

  s = "<table cellpadding=0 cellspacing=0>";

  // write header
  s = s + "<tr>";
  s = s + "<th>gene</th><th>description</th>";
  // FIXME rotate header names
  for(j=0; j<a.arrayName.length; j++) {
    s = s + "<th><small>" + a.arrayName[j].slice(0,2) + "</small></th>";
  }
  s = s + "</tr>";

  // write selected rows of data
  for(i=0; i<30; i++) {
    s = s + "<tr><td>" + a.geneName[i] + "</td><td>" + a.descr[i] + "</td>";

    for(j=0; j<a.x.getSize().width; j++)
      s = s + "<td bgcolor=\"" +
        colorMapping(a.x.getValueAt(i,j)) +
        "\"> </td>";

    s = s + "</tr>";
  }

  s = s + "</table>";

  return s;
}




function showClusters() {


  a = readTSV("cluster.cdt");
  var x1 = standardizeRows(a.x);

  document.getElementById("clusterTable").innerHTML
    = "here are <b>clusters</b> "  +
      showClustersAsTable(a, 0);


}



