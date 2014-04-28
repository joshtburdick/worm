// Utilities for clustering (and drawing clusters.)

goog.require('goog.array');
goog.require('goog.color');
goog.require('goog.math.Matrix');


// Reads in a table written in CDT format.
function readCDT(url) {
  var xmlHttp = null;

  // get entire table, and split by line
  xmlHttp = new XMLHttpRequest();
  // XXX synchronous requests seem to be deprecated; possibly
  // I should be using, e.g., jQuery here?
  xmlHttp.open("get", url, false);
  xmlHttp.overrideMimeType("text/plain")
// FIXME add compression?
//  xmlHttp.setRequestHeaders('accept-encoding', 'gzip');
  xmlHttp.send("");

  var s = xmlHttp.responseText;

  var line = s.split("\n");

  var m = line.length - 3;
  var arrayName = line[0].split("\t").slice(4);
  var geneName = new Array(m);
  var descr = new Array(m);
  var x = new goog.math.Matrix(geneName.length, arrayName.length);

  // loop through the lines, filling in gene info and expression data
  for(i=0; i<m; i++) {
    var f = line[i+3].split("\t");
    geneName[i] = f[1];
    descr[i] = f[2];

    for(j=0; j<arrayName.length; j++)
      x.setValueAt(i, j, parseFloat(f[j+4]));
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
      var a = x.getValueAt(i, j);
      s1 = s1 + a;
      s2 = s2 + a * a;
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

/* Calculates correlation of each gene with a given gene.
  Args:
    selectedIndex - index of gene to center from
  Returns: list with i and foo. */
function calcCorrelations(selectedIndex) {

  var xSelected = new goog.math.Matrix( xStandardized.getSize().width, 1 );
  for(j=0; j<xSelected.getSize().height; j++)
    xSelected.setValueAt(j, 0, xStandardized.getValueAt(selectedIndex, j));

  r = xStandardized.multiply(xSelected);
  r1 = new Array(r.getSize().height);
  var p1 = new Array(xStandardized.getSize().height);
  var n1 = xStandardized.getSize().width - 1;
  for(i=0; i<p1.length; i++) {
    r1[i] = r.getValueAt(i, 0) / n1;
    p1[i] = i;
  }

// alert(r1[0]);
  goog.array.sort(p1, function(a,b) { return r1[b] - r1[a]; } );
  return { r : r1, p : p1 };
}


// Utility to convert a number to a color.
function exprToColor(range, loColor, zeroColor, hiColor) {
  return function (x) {
    if (x >= 0)
      return goog.color.rgbArrayToHex(
        goog.color.blend(hiColor, zeroColor, x / range));
    else
      return goog.color.rgbArrayToHex(
        goog.color.blend(loColor, zeroColor, - x / range));
  };
}

