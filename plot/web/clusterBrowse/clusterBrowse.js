
// goog.require('goog.fs');
goog.require('goog.dom');
goog.require('goog.graphics');
goog.require('goog.math.Matrix');


// Reads in a table written in CDT format.
function readCDT(url) {
  var xmlHttp = null;

  // get entire table, and split by line
  xmlHttp = new XMLHttpRequest();
  xmlHttp.open( "GET", url, false );
// FIXME add compression?
//  xmlHttp.setRequestHeaders('accept-encoding', 'gzip');
  xmlHttp.send( null );
  s = xmlHttp.responseText;
  a = s.split("\n");

  var n = a.length - 3;

  var arrayName = a[0].slice(4);

  var geneName = new Array(n);
  var descr = new Array(n);
  var x = new goog.math.Matrix(n, arrayName.length);

  // loop through the lines, filling in gene info and expression data
  var line;
  for(i=0; i<n; i++) {
    var line = a[i+3].split("\t");
    geneName[i] = line[1];
    descr[i] = line[2];

    for(j=0; j<arrayName.length; j++)
      x.setValueAt(i, j, line[j+3]);
  }

  return { arrayName:arrayName,
    geneName:geneName, descr:descr, x:x };
}

// Computes mean and s.d. of the rows of a matrix.
function rowStats(x) {
  var m = x.getSize().height;
  var n = x.getSize().width;

  var mean = new goog.math.Matrix(m, 1);
  var sd = new goog.math.Matrix(m, 1);

  // loop through the rows
  for(i=0; i<m; i++) {
    var s1 = 0;
    var s2 = 0;

    // add up stats about each entry
    for(j=0; j<n; j++) {
      x1 = x.getValueAt(i, j);
      s1 = s1 + x1;
      s2 = s2 + x1 * x1;
    }

    var m1 = s1 / n;
    mean.setValueAt(i, 0, m1);
    sd.setValueAt(i, 0,
      Math.sqrt( (n/(n-1)) * (s2/n - m1*m1) ) );
  }

  return { mean:mean, sd:sd };
}

// Standardizes the rows of a matrix.
function standardizeRows(x) {
  var x1 = new goog.math.Matrix( x.getSize() );
  var stats = rowStats(x);

  for(i = 0; i<x.getSize().height; i++) {
    var m = stats.mean.getValueAt(i, 0);
    var s = stats.sd.getValueAt(i, 0);
    for(j = 0; j<x.getSize().width; j++)
      x1.setValueAt(i, j, (x.getValueAt(i, j) - m) / s);
  }

  return x1;
}


function showClusters() {


  a = readCDT("cluster.cdt");
  var x1 = standardizeRows(a.x);

  document.getElementById("clusterTable").innerHTML
    = "here are <b>clusters</b>";



}




