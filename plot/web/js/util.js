// Various utility functions.

// Reads a .tsv file, as an array of arrays.
function readTSV(url) {

  var xmlHttp = null;

  xmlHttp = new XMLHttpRequest();
  // XXX synchronous requests seem to be deprecated;
  // possibly I should be using, e.g., jQuery here?
  xmlHttp.open("get", url, false);
  xmlHttp.overrideMimeType("text/plain")
  xmlHttp.send("");

  var s = xmlHttp.responseText;

  var a = s.split("\n");

  for(i=0; i<a.length; i++) {
    a[i] = a[i].split("\t");
  }

  return a;
}

// Utility to color rows alternating colors.
function alternatingFill(i) {
  if (i % 2 == 0)
    return new goog.graphics.SolidFill('#f8f8f8');
  else
    return new goog.graphics.SolidFill('#e8e8e8');
}





